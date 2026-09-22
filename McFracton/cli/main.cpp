// mcf_run - headless entry point.
//
// For now it only provides the regression mode, which pins the simulation to fixed seeds and
// prints every observable at full precision. Comparing that output against the golden files in
// McFracton/tests/golden is what lets the architecture cleanup claim "the physics did not change".
//
//   mcf_run --regression                      print the report to stdout
//   mcf_run --regression --compare <dir>      compare against <dir>/<system>.txt, exit 1 on any diff
//   mcf_run --regression --update <dir>       (re)write the golden files in <dir>
//
// There is also a statistical mode, which is what a change to the random stream has to be checked
// against, since such a change moves every digit while leaving the physics alone:
//
//   mcf_run --stats [--seeds N]               mean and standard error per observable over N chains
//
// The batch-run mode (config file, temperature schedule, Slurm) follows in Phase 4 of the plan.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "AbelianGaugeCube.h"
#include "AbelianGaugeSquare.h"
#include "McMachine.h"
#include "Spiderweb.h"
#include "System.h"
#include "XYSquare.h"

namespace {

// Fixed so a report only depends on the code, never on the machine or the clock. One seed now
// covers the whole run, initial configuration included, because McMachine owns the generator.
constexpr unsigned int kMachineSeed = 6789u;

// A system the reports run over. The factory lets both modes build a fresh system per chain,
// which matters now that the initial configuration is drawn from McMachine's generator.
struct SystemSpec {
	std::string name;
	std::function<std::unique_ptr<System>()> make;
	bool use_overrelaxation = false;
};

std::vector<SystemSpec> Systems()
{
	// Small lattices: the point is reproducibility, not physics quality. QXYSquare is out of scope.
	std::vector<SystemSpec> specs;
	specs.push_back({ "XYSquare", [] { return std::make_unique<XYSquare>(8); }, false });
	specs.push_back({ "AbelianGaugeSquare", [] { return std::make_unique<AbelianGaugeSquare>(6, 6); }, true });
	specs.push_back({ "AbelianGaugeCube", [] { return std::make_unique<AbelianGaugeCube>(4, 4); }, true });
	specs.push_back({ "Spiderweb", [] { return std::make_unique<Spiderweb>(6, 6, 1.0); }, false });
	return specs;
}

const double kTemperatures[] = { 2.0, 1.0, 0.5 };

McMachine::NumericalParams BaseParams(bool use_overrelaxation)
{
	McMachine::NumericalParams params;
	params.delta = 0.1;
	params.updates_per_sweep = 500;
	params.updates_per_overrelaxation = 100;
	params.overrelax = use_overrelaxation;
	return params;
}

std::string RunCase(const SystemSpec& spec)
{
	McMachine::NumericalParams params = BaseParams(spec.use_overrelaxation);
	std::unique_ptr<System> system = spec.make();
	McMachine machine(params, *system, "regression_unused_log.txt", kMachineSeed);

	std::ostringstream out;
	out << std::setprecision(17);

	for (double temperature : kTemperatures)
	{
		for (int sweep = 0; sweep < 20; sweep++)
		{
			machine.Sweep(params.updates_per_sweep, temperature);
			if (spec.use_overrelaxation)
				machine.Overrelax(params.updates_per_overrelaxation);
		}

		// Name-driven, so a system that adds or drops an observable shows up here on its own.
		const std::vector<std::string> names = system->observableNames();
		const std::vector<double> values = system->measure(temperature);
		for (size_t i = 0; i < names.size(); i++)
			out << spec.name << "\tT=" << temperature << '\t' << names[i] << '\t' << values[i] << '\n';
	}

	return out.str();
}

// --- statistical mode -------------------------------------------------------------------------
//
// One chain per seed, cooled through the same temperatures as the regression report. Within a
// chain the samples are correlated, so the error bar is taken across chains: each chain gives one
// mean per observable, and the reported figure is the mean and standard error of those. Two builds
// whose only difference is the random stream must agree here within a few of those error bars.

constexpr int kStatsThermSweeps = 100;
constexpr int kStatsMeasurements = 40;
constexpr int kStatsSweepsPerMeasurement = 5;

std::string RunStats(int n_seeds)
{
	std::ostringstream out;
	out << std::setprecision(10);

	for (const SystemSpec& spec : Systems())
	{
		// [temperature][observable][chain]
		std::vector<std::vector<std::vector<double>>> chain_means(std::size(kTemperatures));
		std::vector<std::string> names;

		for (int chain = 0; chain < n_seeds; chain++)
		{
			McMachine::NumericalParams params = BaseParams(spec.use_overrelaxation);
			std::unique_ptr<System> system = spec.make();
			McMachine machine(params, *system, "stats_unused_log.txt", kMachineSeed + chain);
			names = system->observableNames();

			for (size_t t = 0; t < std::size(kTemperatures); t++)
			{
				const double temperature = kTemperatures[t];

				// Thermalize with the step size still adapting, as StartSimulation does.
				for (int sweep = 0; sweep < kStatsThermSweeps; sweep++)
				{
					machine.Sweep(params.updates_per_sweep, temperature);
					if (spec.use_overrelaxation)
						machine.Overrelax(params.updates_per_overrelaxation);
				}

				// Then measure with delta frozen, again as the real measurement phase does.
				std::vector<double> sums(names.size(), 0.0);
				for (int m = 0; m < kStatsMeasurements; m++)
				{
					for (int sweep = 0; sweep < kStatsSweepsPerMeasurement; sweep++)
						machine.Sweep(params.updates_per_sweep, temperature, false);
					if (spec.use_overrelaxation)
						machine.Overrelax(params.updates_per_overrelaxation);

					const std::vector<double> values = system->measure(temperature);
					for (size_t i = 0; i < sums.size(); i++)
						sums[i] += values[i];
				}

				if (chain_means[t].empty())
					chain_means[t].resize(names.size());
				for (size_t i = 0; i < sums.size(); i++)
					chain_means[t][i].push_back(sums[i] / kStatsMeasurements);
			}
		}

		for (size_t t = 0; t < std::size(kTemperatures); t++)
		{
			for (size_t i = 0; i < names.size(); i++)
			{
				const std::vector<double>& values = chain_means[t][i];
				const int n = (int)values.size();

				double mean = 0.0;
				for (double v : values)
					mean += v;
				mean /= n;

				double error = 0.0;
				if (n > 1)
				{
					double sum_sqr = 0.0;
					for (double v : values)
						sum_sqr += (v - mean) * (v - mean);
					error = std::sqrt(sum_sqr / (double(n) * double(n - 1)));
				}

				out << spec.name << "\tT=" << kTemperatures[t] << '\t' << names[i]
					<< '\t' << mean << '\t' << error << '\t' << n << '\n';
			}
		}

		std::cerr << "stats: " << spec.name << " done\n";
	}

	return out.str();
}

std::string ReadFile(const std::filesystem::path& path)
{
	std::ifstream in(path, std::ios::binary);
	if (!in)
		return {};
	std::ostringstream buffer;
	buffer << in.rdbuf();
	return buffer.str();
}

// Reports the first differing line, which is the one thing you want to see when a refactor slips.
bool ReportDifference(const std::string& name, const std::string& expected, const std::string& actual)
{
	if (expected == actual)
		return false;

	std::istringstream e(expected), a(actual);
	std::string line_e, line_a;
	int line_number = 1;
	for (;;)
	{
		const bool got_e = static_cast<bool>(std::getline(e, line_e));
		const bool got_a = static_cast<bool>(std::getline(a, line_a));
		if (!got_e && !got_a)
			break;
		if (!got_e)
			line_e = "<end of file>";
		if (!got_a)
			line_a = "<end of file>";

		if (line_e != line_a)
		{
			std::cout << "  first difference at line " << line_number << ":\n"
				<< "    golden: " << line_e << "\n"
				<< "    actual: " << line_a << "\n";
			break;
		}
		line_number++;
	}
	return true;
}

int RunRegression(const std::filesystem::path& compare_dir, const std::filesystem::path& update_dir)
{
	int failures = 0;

	for (const SystemSpec& c : Systems())
	{
		const std::string report = RunCase(c);

		if (!update_dir.empty())
		{
			std::filesystem::create_directories(update_dir);
			const std::filesystem::path out_path = update_dir / (c.name + ".txt");
			std::ofstream out(out_path, std::ios::binary);
			out << report;
			std::cout << "wrote " << out_path.string() << '\n';
		}
		else if (!compare_dir.empty())
		{
			const std::filesystem::path golden_path = compare_dir / (c.name + ".txt");
			const std::string golden = ReadFile(golden_path);
			if (golden.empty())
			{
				std::cout << "MISSING " << c.name << " (" << golden_path.string() << ")\n";
				failures++;
			}
			else if (ReportDifference(c.name, golden, report))
			{
				std::cout << "FAIL    " << c.name << '\n';
				failures++;
			}
			else
			{
				std::cout << "ok      " << c.name << '\n';
			}
		}
		else
		{
			std::cout << report;
		}
	}

	if (!compare_dir.empty())
		std::cout << (failures == 0 ? "regression: all systems identical\n"
			: "regression: " + std::to_string(failures) + " system(s) differ\n");

	return failures == 0 ? 0 : 1;
}

} // namespace

int main(int argc, char** argv)
{
	std::vector<std::string> args(argv + 1, argv + argc);

	const char* usage =
		"usage: mcf_run --regression [--compare <dir> | --update <dir>]\n"
		"       mcf_run --stats [--seeds <n>]\n";

	bool regression = false;
	bool stats = false;
	int n_seeds = 10;
	std::filesystem::path compare_dir;
	std::filesystem::path update_dir;

	for (size_t i = 0; i < args.size(); i++)
	{
		const std::string& a = args[i];
		if (a == "--regression")
			regression = true;
		else if (a == "--stats")
			stats = true;
		else if (a == "--seeds" && i + 1 < args.size())
			n_seeds = std::stoi(args[++i]);
		else if (a == "--compare" && i + 1 < args.size())
			compare_dir = args[++i];
		else if (a == "--update" && i + 1 < args.size())
			update_dir = args[++i];
		else
		{
			std::cerr << "unknown argument: " << a << "\n" << usage;
			return 2;
		}
	}

	// Exactly one mode, so a bare run prints the usage rather than silently doing the wrong thing.
	if (regression == stats)
	{
		std::cerr << usage;
		return 2;
	}

	if (stats)
	{
		if (n_seeds < 2)
		{
			std::cerr << "--seeds needs at least 2, so the error bar has something to average\n";
			return 2;
		}
		std::cout << RunStats(n_seeds);
		return 0;
	}

	return RunRegression(compare_dir, update_dir);
}
