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
// The batch-run mode (config file, temperature schedule, Slurm) follows in Phase 4 of the plan.

#include <filesystem>
#include <fstream>
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

// Fixed everywhere so a report only depends on the code, never on the machine or the clock.
constexpr unsigned int kSystemSeed = 12345u;
constexpr unsigned int kMachineSeed = 6789u;

struct RegressionCase {
	std::string name;
	std::unique_ptr<System> system;
	bool use_overrelaxation = false;
};

std::vector<RegressionCase> MakeCases()
{
	// Small lattices: the point is reproducibility, not physics quality. QXYSquare is out of scope.
	std::vector<RegressionCase> cases;
	cases.push_back({ "XYSquare", std::make_unique<XYSquare>(8), false });
	cases.push_back({ "AbelianGaugeSquare", std::make_unique<AbelianGaugeSquare>(6, 6, kSystemSeed), true });
	cases.push_back({ "AbelianGaugeCube", std::make_unique<AbelianGaugeCube>(4, 4, kSystemSeed), true });
	cases.push_back({ "Spiderweb", std::make_unique<Spiderweb>(6, 6, 1.0, kSystemSeed), false });
	return cases;
}

std::string RunCase(RegressionCase& c)
{
	McMachine::NumericalParams params;
	params.delta = 0.1;
	params.updates_per_sweep = 500;
	params.updates_per_overrelaxation = 100;
	params.overrelax = c.use_overrelaxation;

	McMachine machine(params, *c.system, "regression_unused_log.txt", kMachineSeed);

	std::ostringstream out;
	out << std::setprecision(17);

	const double temperatures[] = { 2.0, 1.0, 0.5 };
	for (double temperature : temperatures)
	{
		for (int sweep = 0; sweep < 20; sweep++)
		{
			machine.Sweep(params.updates_per_sweep, temperature);
			if (c.use_overrelaxation)
				machine.Overrelax(params.updates_per_overrelaxation);
		}

		const System::Observables o = c.system->Measure(temperature);
		out << c.name << "\tT=" << temperature << "\tenergy\t" << o.energy << '\n';
		out << c.name << "\tT=" << temperature << "\thelicity_modulus\t" << o.helicity_modulus << '\n';
		out << c.name << "\tT=" << temperature << "\tpolyakov_loop\t" << o.polyakov_loop << '\n';
		out << c.name << "\tT=" << temperature << "\tflux_cos\t" << o.flux_cos << '\n';
		out << c.name << "\tT=" << temperature << "\tn_defects_a\t" << o.n_defects_a << '\n';
		out << c.name << "\tT=" << temperature << "\tn_defects_b\t" << o.n_defects_b << '\n';
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
	std::vector<RegressionCase> cases = MakeCases();
	int failures = 0;

	for (RegressionCase& c : cases)
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

	bool regression = false;
	std::filesystem::path compare_dir;
	std::filesystem::path update_dir;

	for (size_t i = 0; i < args.size(); i++)
	{
		const std::string& a = args[i];
		if (a == "--regression")
			regression = true;
		else if (a == "--compare" && i + 1 < args.size())
			compare_dir = args[++i];
		else if (a == "--update" && i + 1 < args.size())
			update_dir = args[++i];
		else
		{
			std::cerr << "unknown argument: " << a << "\n"
				<< "usage: mcf_run --regression [--compare <dir> | --update <dir>]\n";
			return 2;
		}
	}

	if (!regression)
	{
		std::cerr << "usage: mcf_run --regression [--compare <dir> | --update <dir>]\n";
		return 2;
	}

	return RunRegression(compare_dir, update_dir);
}
