#include <SFML/Graphics.hpp>
#include <imgui.h>
#include <imgui-SFML.h>

#include <iostream>
#include <algorithm>

#include "Canvas.h"
#include "XYSquare.h"
#include "QXYSquare.h"
#include "AbelianGaugeSquare.h"
#include "Timer.h"
#include "McMachine.h"

int main() {

	using namespace sf;

	const int tau_layers = 15;
	AbelianGaugeSquare cubicLattice(15, tau_layers);
	McMachine::NumericalParams params;
	params.t_max = 100.0;
	params.measure_sweeps = 500;
	params.t_min = 0.01;
	params.max_therm_sweeps = 2000;
	params.overrelax = false;
	McMachine machine(params, cubicLattice, "data.txt");
	
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
	canvas.Initialize(cubicLattice, 15.0);

	bool pause = true;
	bool draw_fluxes = false;
	int field_direction = 0;
	bool plot_energies = false;

	sf::Clock clock;
	Int32 elapsedTime = 0;

	float temperature = 0.1;
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

				if (e.key.code == Keyboard::Space)
				{
					pause = !pause;
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

		ImGui::Begin("Monte Carlo Simulation");
		ImGui::SliderFloat("Temperature", &temperature, 0.0001f, 10.0f, "%.5f");
		//ImGui::SliderFloat("K_s", &squareLattice.K_s, 0.1f, 10.0, "%.3f");
		//ImGui::SliderFloat("K_t", &squareLattice.K_t, 0.1f, 10.0, "%.3f");
		ImGui::SliderInt("Field Direction", &field_direction, 0, 2);
		ImGui::Checkbox("Draw Fluxes", &draw_fluxes);
		ImGui::Checkbox("Plot Energy", &plot_energies);
		ImGui::Checkbox("Pause", &pause);
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

		//if(!draw_s_vortices && !draw_t_vortices)
		//	canvas.Draw(squareLattice, layer);
		//if (draw_s_vortices)
		//	canvas.Draw(squareLattice, squareLattice.getSpacialVortices(), layer);
		//else if (draw_t_vortices)
		//	canvas.Draw(squareLattice, squareLattice.getTemporalVortices(), layer);

		//int n_monopoles = 0;
		//const auto monopoles = cubicLattice.getMonopoles();
		//std::for_each(monopoles.begin(), monopoles.end(), [&n_monopoles](int m) { n_monopoles += std::abs(m); });
		//std::cout << n_monopoles << std::endl;

		if (draw_fluxes)
			canvas.DrawFluxes(cubicLattice, layer);
		else
			canvas.Draw(cubicLattice, field_direction, layer);

		ImGui::SFML::Render(window);
		window.display();

		if (!pause) {
			machine.Sweep(2000, temperature);
			//machine.Overrelax(5000);
		}

		double current_energy = cubicLattice.getEnergy();
		energies.Push((float)current_energy);
	}

	return 0;
}