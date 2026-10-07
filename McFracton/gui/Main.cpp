// mcf_gui - watch a simulation run.
//
// The window is driven entirely by the registry: pick a system, set its parameters, press Rebuild.
// Nothing here names a system class, so a system added to core/include/SystemRegistry.h shows up in
// the picker on its own, and LatticeView draws whatever channels it declares.
//
// Annealing runs belong to mcf_run; this is for looking at a system at one temperature.
//
//   mcf_gui [system]      start on that registered system instead of Spiderweb_corrected

#include <SFML/Graphics.hpp>
#include <imgui.h>
#include <imgui-SFML.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "BufferedArray.h"
#include "LatticeView.h"
#include "McMachine.h"
#include "System.h"
#include "SystemRegistry.h"

#ifdef _WIN32
#include <windows.h>
#endif

int main(int argc, char** argv)
{
	// Without this, Windows hands a scaled display (150% here) a coordinate space smaller than the
	// screen while SFML still reports the screen's true size, so a window sized from the latter
	// overflows by exactly the scale factor and the right of the lattice falls off the edge.
	// Being DPI aware makes the two agree, and renders at native resolution rather than upscaling.
#ifdef _WIN32
	SetProcessDPIAware();
#endif

	const std::vector<mcf::SystemEntry>& registry = mcf::systemRegistry();

	// Starts on Spiderweb_corrected unless the first argument names another system: mcf_gui XYSquare
	const std::string start_system = argc > 1 ? argv[1] : "Spiderweb_corrected";
	int selected = 0;
	for (size_t i = 0; i < registry.size(); i++)
		if (registry[i].name == start_system)
			selected = (int)i;
	std::vector<double> parameters;
	for (const mcf::Parameter& p : registry[selected].parameters)
		parameters.push_back(p.default_value);

	McMachine::NumericalParams params;

	std::unique_ptr<System> system = registry[selected].make(parameters);
	std::unique_ptr<McMachine> machine = std::make_unique<McMachine>(params, *system);

	// 1900x1200 where the desktop allows it, and never larger than the desktop, since a window
	// bigger than the screen silently hides the right hand side of the lattice.
	const sf::VideoMode desktop = sf::VideoMode::getDesktopMode();
	const unsigned int width = std::min(1900u, (unsigned int)(desktop.width * 0.95));
	const unsigned int height = std::min(1200u, (unsigned int)(desktop.height * 0.90));

	std::cout << "desktop " << desktop.width << "x" << desktop.height
		<< ", window " << width << "x" << height << std::endl;

	sf::RenderWindow window(sf::VideoMode(width, height), "McFracton");
	window.setVerticalSyncEnabled(true);
	if (!ImGui::SFML::Init(window))
		return -1;

	ImGuiIO& io = ImGui::GetIO();
	io.Fonts->Clear();
	io.Fonts->AddFontFromFileTTF("DMSans-VariableFont_opsz,wght.ttf", 22.0f);
	ImGui::SFML::UpdateFontTexture();

	ImGui::GetStyle().ScaleAllSizes(1.5f);

	const float panel_width = 640.0f;
	LatticeView view(window);
	view.setSystem(*system);
	view.setArea(panel_width, 0.0f, (float)window.getSize().x - panel_width, (float)window.getSize().y);

	bool pause = true;
	std::string build_error;
	float temperature = 1.0f;
	int updates_per_frame = 2000;
	bool plot_energy = true;
	BufferedArray energies(200);

	sf::Clock clock;
	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			ImGui::SFML::ProcessEvent(window, event);

			if (event.type == sf::Event::Closed)
				window.close();

			if (event.type == sf::Event::Resized)
				view.setArea(panel_width, 0.0f, (float)event.size.width - panel_width,
					(float)event.size.height);

			// The keys of the old canvases: the arrows step the first sliced axis, Enter the
			// second, Space pauses. Up and down step the channel, and with Ctrl held the overlay.
			// ImGui gets first refusal so typing in a field still works.
			if (event.type == sf::Event::KeyPressed && !ImGui::GetIO().WantCaptureKeyboard)
			{
				switch (event.key.code)
				{
				case sf::Keyboard::Right: view.stepSlice(*system, 0, +1); break;
				case sf::Keyboard::Left:  view.stepSlice(*system, 0, -1); break;
				case sf::Keyboard::Enter: view.stepSlice(*system, 1, +1); break;
				case sf::Keyboard::Space: pause = !pause; break;
				// Down moves down the list in the combo, which is the next channel.
				case sf::Keyboard::Down:
					if (event.key.control) view.stepOverlay(*system, +1);
					else                   view.stepChannel(*system, +1);
					break;
				case sf::Keyboard::Up:
					if (event.key.control) view.stepOverlay(*system, -1);
					else                   view.stepChannel(*system, -1);
					break;
				default: break;
				}
			}
		}

		ImGui::SFML::Update(window, clock.restart());

		// The panel is pinned to the strip the lattice view leaves free, so the two never overlap
		// and the labels always have room.
		ImGui::SetNextWindowPos({ 0.0f, 0.0f });
		ImGui::SetNextWindowSize({ panel_width, (float)window.getSize().y });
		ImGui::Begin("McFracton", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize
			| ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar);
		// Widgets take the left half, leaving the right half for the label.
		ImGui::PushItemWidth(panel_width * 0.42f);

		// --- the system picker, straight out of the registry ---
		std::vector<const char*> names;
		for (const mcf::SystemEntry& entry : registry)
			names.push_back(entry.name.c_str());
		if (ImGui::Combo("System", &selected, names.data(), (int)names.size()))
		{
			parameters.clear();
			for (const mcf::Parameter& p : registry[selected].parameters)
				parameters.push_back(p.default_value);
		}

		for (size_t i = 0; i < registry[selected].parameters.size(); i++)
		{
			const mcf::Parameter& p = registry[selected].parameters[i];
			if (p.is_integer)
			{
				int value = (int)parameters[i];
				if (ImGui::InputInt(p.name.c_str(), &value))
					parameters[i] = std::max(1, value);
			}
			else
			{
				float value = (float)parameters[i];
				if (ImGui::InputFloat(p.name.c_str(), &value, 0.0f, 0.0f, "%.4f"))
					parameters[i] = value;
			}
		}

		if (ImGui::Button("Rebuild"))
		{
			// A system may refuse its parameters; then the one on screen simply stays.
			try
			{
				std::unique_ptr<System> rebuilt = registry[selected].make(parameters);
				// The machine holds a reference to the system, so it goes first.
				machine.reset();
				system = std::move(rebuilt);
				machine = std::make_unique<McMachine>(params, *system);
				view.setSystem(*system);
				energies = BufferedArray(200);
				pause = true;
				build_error.clear();
			}
			catch (const std::invalid_argument& e)
			{
				build_error = e.what();
			}
		}
		if (!build_error.empty())
			ImGui::TextWrapped("%s", build_error.c_str());

		ImGui::Separator();

		// --- running ---
		ImGui::SliderFloat("Temperature", &temperature, 0.01f, 10.0f, "%.5f");
		ImGui::SliderInt("Updates / frame", &updates_per_frame, 100, 50000);
		ImGui::Checkbox("Pause", &pause);
		ImGui::SameLine();
		ImGui::TextDisabled("(space)");

		ImGui::Separator();

		// --- what is on screen ---
		view.drawControls(*system);

		ImGui::Separator();

		ImGui::Checkbox("Plot energy", &plot_energy);
		if (plot_energy)
		{
			ImGui::PlotLines("Energy", energies.get_data().data(), energies.get_size(),
				energies.get_offset(), nullptr, energies.get_min(), energies.get_max(),
				ImVec2(0, 150));
		}
		ImGui::Text("%d variables", system->numVariables());

		ImGui::PopItemWidth();
		ImGui::End();

		if (!pause)
			machine->Sweep(updates_per_frame, temperature);

		energies.Push((float)system->getEnergy());

		window.clear();
		view.draw(*system);
		ImGui::SFML::Render(window);
		window.display();
	}

	ImGui::SFML::Shutdown();
	return 0;
}
