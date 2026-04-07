#include <SFML/Graphics.hpp>
#include <imgui.h>
#include <imgui-SFML.h>

#include <iostream>

#include "Canvas.h"
#include "XYSquare.h"
#include "Timer.h"
#include "McMachine.h"

int main() {

	using namespace sf;

	XYSquare squareLattice(80);
	McMachine::NumericalParams params;
	McMachine machine(params, squareLattice, "data.txt");

	machine.StartSimulation();

	RenderWindow window(VideoMode(1900, 1200), "Simulation");
	window.setVerticalSyncEnabled(true);
	if (!ImGui::SFML::Init(window))
		return -1;

	ImGuiIO& io = ImGui::GetIO();
	io.Fonts->Clear();
	io.Fonts->AddFontFromFileTTF("DMSans-VariableFont_opsz,wght.ttf", 22.0f);
	ImGui::SFML::UpdateFontTexture();

	ImGuiStyle& style = ImGui::GetStyle();
	style.ScaleAllSizes(1.5f);

	Canvas canvas(window, { 200.0f, 0.0f });
	canvas.Initialize(squareLattice, 12.0);

	bool draw_vortices = false;
	bool plot_energies = false;

	sf::Clock clock;
	Int32 elapsedTime = 0;

	float temperature = 10.0;
	BufferedArray energies(100);
	while (window.isOpen())
	{
		Event e;
		while (window.pollEvent(e))
		{
			ImGui::SFML::ProcessEvent(window, e);

			if (e.type == Event::Closed)
			{
				window.close();
			}
		}

		sf::Time dt = clock.restart();
		elapsedTime += dt.asMilliseconds();

		ImGui::SFML::Update(window, dt);

		ImGui::Begin("Hello, world!");
		ImGui::SliderFloat("Temperature", &temperature, 0.01f, 10.0f, "%.3f");
		ImGui::Checkbox("Vortices", &draw_vortices);
		ImGui::Checkbox("Plot Energy", &plot_energies);
		if(plot_energies) {
			ImGui::PlotLines(
				"Energy",
				energies.get_data().data(),
				energies.get_size(),
				energies.get_offset(),
				nullptr,
				energies.get_min(), energies.get_max(),
				ImVec2(0, 150)
			);
		}
		ImGui::End();

		window.clear();

		if (draw_vortices)
			canvas.Draw(squareLattice, squareLattice.getVortices());
		else
			canvas.Draw(squareLattice);

		ImGui::SFML::Render(window);
		window.display();

		machine.Sweep(5000, (double)temperature);

		double current_energy = squareLattice.getEnergy();
		energies.Push((float)current_energy);
	}

	return 0;
}