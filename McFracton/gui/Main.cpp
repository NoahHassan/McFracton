#include <SFML/Graphics.hpp>
#include <imgui.h>
#include <imgui-SFML.h>

#include <iostream>

#include "Canvas.h"
#include "XYSquare.h"
#include "QXYSquare.h"
#include "Timer.h"
#include "McMachine.h"

int main() {

	using namespace sf;

	const int tau_layers = 30;
	QXYSquare squareLattice(30, tau_layers, 1.0, 1.0);
	McMachine::NumericalParams params;
	params.overrelax = false;
	McMachine machine(params, squareLattice, "data.txt");
	
	//machine.StartSimulation();

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

	bool draw_s_vortices = false;
	bool draw_t_vortices = false;
	bool plot_energies = false;

	sf::Clock clock;
	Int32 elapsedTime = 0;

	float temperature = 1.0;
	int layer = 0;
	BufferedArray energies(200);
	while (window.isOpen())
	{
		Event e;
		while (window.pollEvent(e))
		{
			ImGui::SFML::ProcessEvent(window, e);

			if (e.type == Event::KeyPressed)
			{
				if (e.key.code == Keyboard::Enter)
				{
					layer = (layer + 1) % tau_layers;
				}
			}

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
		ImGui::SliderFloat("K_s", &squareLattice.K_s, 0.1f, 10.0, "%.3f");
		ImGui::SliderFloat("K_t", &squareLattice.K_t, 0.1f, 10.0, "%.3f");
		if (ImGui::Checkbox("Spacial Vortices", &draw_s_vortices))
		{
			draw_t_vortices = false;
		}
		if (ImGui::Checkbox("Temporal Vortices", &draw_t_vortices))
		{
			draw_s_vortices = false;
		}
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

		if(!draw_s_vortices && !draw_t_vortices)
			canvas.Draw(squareLattice, layer);
		if (draw_s_vortices)
			canvas.Draw(squareLattice, squareLattice.getSpacialVortices(), layer);
		else if (draw_t_vortices)
			canvas.Draw(squareLattice, squareLattice.getTemporalVortices(), layer);

		ImGui::SFML::Render(window);
		window.display();

		machine.Sweep(2000, temperature);
		//machine.Overrelax(5000);

		double current_energy = squareLattice.getEnergy();
		energies.Push((float)current_energy);
	}

	return 0;
}