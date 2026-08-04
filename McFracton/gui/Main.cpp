#include <SFML/Graphics.hpp>
#include <imgui.h>
#include <imgui-SFML.h>

#include <iostream>
#include <algorithm>

#include "Canvas.h"
#include "HyperCanvas.h"
#include "XYSquare.h"
#include "QXYSquare.h"
#include "AbelianGaugeSquare.h"
#include "AbelianGaugeCube.h"
#include "Timer.h"
#include "McMachine.h"

int main() {

	using namespace sf;

	const int space_layers = 6;
	const int tau_layers = 6;
	AbelianGaugeCube hypercubicLattice(space_layers, tau_layers);
	McMachine::NumericalParams params;
	params.t_max = 100.0;
	params.t_min = 0.01;
	params.max_therm_sweeps = 4000;
	params.n_measurements = 10;
	params.max_measure_sweeps = 1000;
	params.overrelax = true;
	params.updates_per_overrelaxation = 1000;
	McMachine machine(params, hypercubicLattice, "3d_abelian_L=6.txt");

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

	HyperCanvas canvas(window, { 200.0f, 0.0f });
	canvas.Initialize(hypercubicLattice, 15.0);

	bool pause = true;
	bool draw_monopoles = false;
	bool draw_fluxes = false;
	int field_direction = 0;
	bool plot_energies = false;

	sf::Clock clock;
	Int32 elapsedTime = 0;

	float temperature = 0.1f;
	int layer = 0;
	int time = 0;
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
					layer = (layer + 1) % space_layers;
				}
				if (e.key.code == Keyboard::Right)
				{
					time = (time + 1) % tau_layers;
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
		ImGui::SliderFloat("Temperature", &temperature, 0.01f, 10.0f, "%.5f");
		//ImGui::SliderFloat("K_s", &squareLattice.K_s, 0.1f, 10.0, "%.3f");
		//ImGui::SliderFloat("K_t", &squareLattice.K_t, 0.1f, 10.0, "%.3f");
		ImGui::SliderInt("Field Direction", &field_direction, 0, 2);
		if (ImGui::Checkbox("Draw Fluxes", &draw_fluxes))
		{
			draw_monopoles = false;
		}
		else if (ImGui::Checkbox("Draw Monopoles", &draw_monopoles))
		{
			draw_fluxes = false;
		}
		ImGui::Checkbox("Plot Energy", &plot_energies);
		ImGui::Checkbox("Pause", &pause);
		ImGui::Checkbox("Overrelax", &params.overrelax);
		if (plot_energies) {
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

		if (draw_monopoles)
			canvas.DrawMonopoles(hypercubicLattice, layer, time);
		else if (draw_fluxes)
			canvas.DrawFluxes(hypercubicLattice, layer, time);
		else
			canvas.Draw(hypercubicLattice, field_direction, layer, time);

		ImGui::SFML::Render(window);
		window.display();

		if (!pause) {
			machine.Sweep(2000, temperature);
			if (params.overrelax)
				machine.Overrelax(500);
		}

		//double current_energy = (cubicLattice.getEnergy() / temperature) / (cubicLattice.nPlaqs);
		double current_energy = (hypercubicLattice.getEnergy());
		energies.Push((float)current_energy);
	}

	return 0;
}