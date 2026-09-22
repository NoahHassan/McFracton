// mcf_run - headless entry point.
//
// An annealing run is configured by key=value settings, which come from a config file, from
// --key=value on the command line, or both; the command line wins, so one config file plus a few
// overrides covers a whole Slurm array. Keys are the system name, that system's own constructor
// parameters (see --list), every McMachine::NumericalParams field, and seed / out / git_hash.
//
//   mcf_run --config run.cfg                  run with the settings in run.cfg
//   mcf_run --config run.cfg --seed=7         ... with seed overridden
//   mcf_run --system=Spiderweb --KU=0.5       run with no config file at all
//   mcf_run --list                            print the registered systems and their parameters
//
// Two verification modes share the binary. The regression mode pins the simulation to a fixed seed
// and prints every observable at full precision; comparing that against the golden files in
// McFracton/tests/golden is what lets the architecture cleanup claim "the physics did not change":
//
//   mcf_run --regression                      print the report to stdout
//   mcf_run --regression --compare <dir>      compare against <dir>/<system>.txt, exit 1 on any diff
//   mcf_run --regression --update <dir>       (re)write the golden files in <dir>
//
// The statistical mode is what a change to the random stream has to be checked against instead,
// since such a change moves every digit while leaving the physics alone:
//
//   mcf_run --stats [--seeds N]               mean and standard error per observable over N chains

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "McMachine.h"
#include "System.h"
#include "SystemRegistry.h"

namespace {

// Fixed so a report only depends on the code, never on the machine or the clock. One seed now
// covers the whole run, initial configuration included, because McMachine owns the generator.
constexpr unsigned int kMachineSeed = 6789u;

// A system the reports run over, named the way the registry names it, with the constructor
// parameters spelled out rather than defaulted so the report never moves when a default does.
// Both modes rebuild the system per chain, which matters now that the initial configuration is
// drawn from McMachine's generator.
struct SystemSpec {
	std::string name;
	std::vector<double> parameters;
	bool use_overrelaxation = false;
};

std::unique_ptr<System> Build(const SystemSpec& spec)
{
	const mcf::SystemEntry* entry = mcf::findSystem(spec.name);
	if (entry == nullptr)
		throw std::runtime_error("no registered system named " + spec.name);
	return entry->make(spec.parameters);
}

std::vector<SystemSpec> Systems()
{
	// Small lattices: the point is reproducibility, not physics quality. QXYSquare is out of scope.
	std::vector<SystemSpec> specs;
	specs.push_back({ "XYSquare", { 8 }, false });
	specs.push_back({ "AbelianGaugeSquare", { 6, 6 }, true });
	specs.push_back({ "AbelianGaugeCube", { 4, 4 }, true });
	specs.push_back({ "Spiderweb", { 6, 6, 1.0 }, false });
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
	std::unique_ptr<System> system = Build(spec);
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
			std::unique_ptr<System> system = Build(spec);
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

// --- configured annealing run -----------------------------------------------------------------
//
// Settings are a flat map of strings. A config file supplies some, the command line overrides
// them, and whatever is left keeps the default from NumericalParams or from the registry.

using Settings = std::vector<std::pair<std::string, std::string>>;

std::string Trim(const std::string& text)
{
	const size_t first = text.find_first_not_of(" \t\r\n");
	if (first == std::string::npos)
		return {};
	const size_t last = text.find_last_not_of(" \t\r\n");
	return text.substr(first, last - first + 1);
}

// Later entries win, so the command line can simply be appended to the file.
const std::string* Lookup(const Settings& settings, const std::string& key)
{
	const std::string* found = nullptr;
	for (const auto& entry : settings)
		if (entry.first == key)
			found = &entry.second;
	return found;
}

// A '#' starts a comment; blank lines are skipped. Anything else must be key = value.
bool ReadConfig(const std::filesystem::path& path, Settings& settings, std::string& error)
{
	std::ifstream in(path);
	if (!in)
	{
		error = "cannot open config file " + path.string();
		return false;
	}

	std::string line;
	int line_number = 0;
	while (std::getline(in, line))
	{
		line_number++;
		const size_t comment = line.find('#');
		if (comment != std::string::npos)
			line = line.substr(0, comment);
		line = Trim(line);
		if (line.empty())
			continue;

		const size_t equals = line.find('=');
		if (equals == std::string::npos)
		{
			error = path.string() + ":" + std::to_string(line_number) + ": expected key = value";
			return false;
		}
		settings.emplace_back(Trim(line.substr(0, equals)), Trim(line.substr(equals + 1)));
	}
	return true;
}

// Every key a run understands, so an unknown one can be reported instead of silently ignored:
// a typo in a Slurm config would otherwise cost a whole array job.
std::vector<std::string> KnownKeys(const mcf::SystemEntry& entry)
{
	std::vector<std::string> keys = { "system", "seed", "out", "git_hash",
		"t_max", "t_min", "t_fac", "delta", "updates_per_sweep", "initial_therm_sweeps",
		"max_therm_sweeps", "max_measure_sweeps", "n_measurements", "overrelax",
		"log_energies", "updates_per_overrelaxation" };
	for (const mcf::Parameter& p : entry.parameters)
		keys.push_back(p.name);
	return keys;
}

bool ParseBool(const std::string& text, bool& out, std::string& error)
{
	if (text == "true" || text == "1" || text == "yes") { out = true;  return true; }
	if (text == "false" || text == "0" || text == "no") { out = false; return true; }
	error = "expected true or false, got " + text;
	return false;
}

// Reads one setting into `target`, leaving it alone when the key is absent.
template <typename T>
bool ReadNumber(const Settings& settings, const std::string& key, T& target, std::string& error)
{
	const std::string* text = Lookup(settings, key);
	if (text == nullptr)
		return true;
	try
	{
		target = (T)std::stod(*text);
	}
	catch (const std::exception&)
	{
		error = key + ": expected a number, got " + *text;
		return false;
	}
	return true;
}

void PrintRegistry()
{
	for (const mcf::SystemEntry& entry : mcf::systemRegistry())
	{
		std::cout << entry.name << '\n';
		for (const mcf::Parameter& p : entry.parameters)
		{
			std::cout << "    " << p.name << " = ";
			if (p.is_integer)
				std::cout << (long long)p.default_value;
			else
				std::cout << p.default_value;
			std::cout << '\n';
		}
	}
}

int RunConfigured(const Settings& settings)
{
	std::string error;

	const std::string* system_name = Lookup(settings, "system");
	if (system_name == nullptr)
	{
		std::cerr << "no system given: set system = <name>, or run --list to see the choices\n";
		return 2;
	}
	const mcf::SystemEntry* entry = mcf::findSystem(*system_name);
	if (entry == nullptr)
	{
		std::cerr << "unknown system " << *system_name << "; run --list to see the choices\n";
		return 2;
	}

	const std::vector<std::string> known = KnownKeys(*entry);
	for (const auto& setting : settings)
	{
		if (std::find(known.begin(), known.end(), setting.first) == known.end())
		{
			std::cerr << "unknown setting " << setting.first << " for system " << entry->name
				<< "; run --list to see that system's parameters\n";
			return 2;
		}
	}

	// System parameters, each defaulted by the registry unless the settings name it.
	std::vector<double> values;
	for (const mcf::Parameter& p : entry->parameters)
	{
		double value = p.default_value;
		if (!ReadNumber(settings, p.name, value, error))
		{
			std::cerr << error << '\n';
			return 2;
		}
		values.push_back(value);
	}

	McMachine::NumericalParams params;
	const bool ok = ReadNumber(settings, "t_max", params.t_max, error)
		&& ReadNumber(settings, "t_min", params.t_min, error)
		&& ReadNumber(settings, "t_fac", params.t_fac, error)
		&& ReadNumber(settings, "delta", params.delta, error)
		&& ReadNumber(settings, "updates_per_sweep", params.updates_per_sweep, error)
		&& ReadNumber(settings, "initial_therm_sweeps", params.initial_therm_sweeps, error)
		&& ReadNumber(settings, "max_therm_sweeps", params.max_therm_sweeps, error)
		&& ReadNumber(settings, "max_measure_sweeps", params.max_measure_sweeps, error)
		&& ReadNumber(settings, "n_measurements", params.n_measurements, error)
		&& ReadNumber(settings, "updates_per_overrelaxation", params.updates_per_overrelaxation, error);
	if (!ok)
	{
		std::cerr << error << '\n';
		return 2;
	}

	if (const std::string* text = Lookup(settings, "overrelax"))
		if (!ParseBool(*text, params.overrelax, error))
		{
			std::cerr << "overrelax: " << error << '\n';
			return 2;
		}
	if (const std::string* text = Lookup(settings, "log_energies"))
		if (!ParseBool(*text, params.log_energies, error))
		{
			std::cerr << "log_energies: " << error << '\n';
			return 2;
		}

	unsigned int seed = 0;   // 0 still means "draw a fresh one and record it".
	if (!ReadNumber(settings, "seed", seed, error))
	{
		std::cerr << error << '\n';
		return 2;
	}

	// Runs land in results/ by default, which .gitignore already covers, rather than scattering
	// logs through the repository root. McMachine asserts on a log it cannot open, so the
	// directory is created first - on a cluster that assert would be a puzzling way to fail.
	const std::string* out = Lookup(settings, "out");
	const std::filesystem::path logfile =
		out != nullptr ? std::filesystem::path(*out)
		               : std::filesystem::path("results") / (entry->name + ".txt");
	if (logfile.has_parent_path() && !logfile.parent_path().empty())
	{
		std::error_code ec;
		std::filesystem::create_directories(logfile.parent_path(), ec);
		if (ec)
		{
			std::cerr << "cannot create " << logfile.parent_path().string() << ": " << ec.message() << '\n';
			return 2;
		}
	}

	std::unique_ptr<System> system = entry->make(values);
	McMachine machine(params, *system, logfile.string(), seed);

	// Provenance McMachine cannot know by itself. The git hash is passed in rather than compiled
	// in, so neither build system needs a generated header; the Slurm script fills it from git.
	machine.addProvenance("system", entry->name);
	for (size_t i = 0; i < entry->parameters.size(); i++)
	{
		std::ostringstream value;
		value << std::setprecision(17) << values[i];
		machine.addProvenance(entry->parameters[i].name, value.str());
	}
	if (const std::string* git_hash = Lookup(settings, "git_hash"))
		if (!git_hash->empty())
			machine.addProvenance("git_hash", *git_hash);

	std::cout << "running " << entry->name << " -> " << logfile.string() << '\n';
	machine.StartSimulation();
	std::cout << "wrote " << logfile.string() << '\n';
	return 0;
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
		"usage: mcf_run --config <file> [--key=value ...]\n"
		"       mcf_run --key=value ...\n"
		"       mcf_run --list\n"
		"       mcf_run --regression [--compare <dir> | --update <dir>]\n"
		"       mcf_run --stats [--seeds <n>]\n";

	bool regression = false;
	bool stats = false;
	bool list = false;
	int n_seeds = 10;
	std::filesystem::path compare_dir;
	std::filesystem::path update_dir;
	std::filesystem::path config_path;
	Settings overrides;

	for (size_t i = 0; i < args.size(); i++)
	{
		const std::string& a = args[i];
		if (a == "--regression")
			regression = true;
		else if (a == "--stats")
			stats = true;
		else if (a == "--list")
			list = true;
		else if (a == "--seeds" && i + 1 < args.size())
			n_seeds = std::stoi(args[++i]);
		else if (a == "--compare" && i + 1 < args.size())
			compare_dir = args[++i];
		else if (a == "--update" && i + 1 < args.size())
			update_dir = args[++i];
		else if (a == "--config" && i + 1 < args.size())
			config_path = args[++i];
		else if (a.rfind("--", 0) == 0 && a.find('=') != std::string::npos)
		{
			// --key=value, the form a Slurm array uses to vary one setting per task.
			const size_t equals = a.find('=');
			overrides.emplace_back(Trim(a.substr(2, equals - 2)), Trim(a.substr(equals + 1)));
		}
		else
		{
			std::cerr << "unknown argument: " << a << "\n" << usage;
			return 2;
		}
	}

	if (list)
	{
		PrintRegistry();
		return 0;
	}

	// The verification modes are self-contained, so they refuse to be mixed with run settings.
	if (regression || stats)
	{
		if (regression && stats)
		{
			std::cerr << usage;
			return 2;
		}
		if (!config_path.empty() || !overrides.empty())
		{
			std::cerr << "--regression and --stats take their own fixed settings; "
				"drop the config and the overrides\n";
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

	Settings settings;
	if (!config_path.empty())
	{
		std::string error;
		if (!ReadConfig(config_path, settings, error))
		{
			std::cerr << error << '\n';
			return 2;
		}
	}
	// Appended last, so a command-line override beats the same key in the file.
	settings.insert(settings.end(), overrides.begin(), overrides.end());

	if (settings.empty())
	{
		std::cerr << usage;
		return 2;
	}

	return RunConfigured(settings);
}
