#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "AbelianGaugeCube.h"
#include "AbelianGaugeSquare.h"
#include "QXYSquare.h"
#include "Spiderweb.h"
#include "System.h"
#include "XYSquare.h"

/// <summary>
/// The one place a system has to announce itself. A name, the parameters its constructor needs,
/// and a factory that turns values for those parameters into a system.
///
/// Everything that has to offer a choice of system reads this list instead of naming the classes:
/// the CLI resolves config keys against it, and the GUI builds its picker from it. Adding a system
/// then means adding one line here, and nothing in cli/ or gui/ has to be touched.
/// </summary>
namespace mcf {

/// <summary>
/// One constructor argument. Values travel as doubles because that is enough for every parameter
/// the systems take; is_integer only says how it should be read and shown.
/// </summary>
struct Parameter {
	std::string name;
	double default_value;
	bool is_integer;
};

struct SystemEntry {
	std::string name;
	std::vector<Parameter> parameters;
	// Values are in the same order as `parameters`, already defaulted by the caller.
	std::function<std::unique_ptr<System>(const std::vector<double>&)> make;
};

inline const std::vector<SystemEntry>& systemRegistry()
{
	static const std::vector<SystemEntry> registry = {
		{
			"XYSquare",
			{ { "size", 16, true } },
			[](const std::vector<double>& p) -> std::unique_ptr<System> {
				return std::make_unique<XYSquare>((int)p[0]);
			}
		},
		{
			"QXYSquare",
			{ { "size", 16, true }, { "Ntau", 16, true }, { "K_s", 1.0, false }, { "K_t", 1.0, false } },
			[](const std::vector<double>& p) -> std::unique_ptr<System> {
				return std::make_unique<QXYSquare>((int)p[0], (int)p[1], (float)p[2], (float)p[3]);
			}
		},
		{
			"AbelianGaugeSquare",
			{ { "linear_size", 8, true }, { "temporal_size", 8, true } },
			[](const std::vector<double>& p) -> std::unique_ptr<System> {
				return std::make_unique<AbelianGaugeSquare>((int)p[0], (int)p[1]);
			}
		},
		{
			"AbelianGaugeCube",
			{ { "linear_size", 6, true }, { "temporal_size", 6, true } },
			[](const std::vector<double>& p) -> std::unique_ptr<System> {
				return std::make_unique<AbelianGaugeCube>((int)p[0], (int)p[1]);
			}
		},
		{
			"Spiderweb",
			{ { "linear_size", 16, true }, { "temporal_size", 16, true }, { "KU", 0.5, false } },
			[](const std::vector<double>& p) -> std::unique_ptr<System> {
				return std::make_unique<Spiderweb>((int)p[0], (int)p[1], p[2]);
			}
		},
	};
	return registry;
}

/// <summary>
/// The entry with this name, or nullptr if there is none.
/// </summary>
inline const SystemEntry* findSystem(const std::string& name)
{
	for (const SystemEntry& entry : systemRegistry())
		if (entry.name == name)
			return &entry;
	return nullptr;
}

} // namespace mcf
