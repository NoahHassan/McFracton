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

	XYSquare squareLattice(80);
	Canvas canvas(window, { 200.0f, 0.0f });
	canvas.Initialize(squareLattice, 12.0);

	McMachine::NumericalParams params(
		100.0, 0.1, 0.1
	);
	McMachine machine(params, squareLattice);

	bool draw_vortices = false;

	sf::Clock clock;
	Int32 elapsedTime = 0;

	float temperature = 10.0;
	while (window.isOpen())
	{
		Event e;
		while (window.pollEvent(e))
		{
			ImGui::SFML::ProcessEvent(window, e);
			if (e.type == Event::KeyPressed)
			{
				if (e.key.code == Keyboard::Key::Enter)
				{
					std::cout << squareLattice.getEnergy() << std::endl;
				}

				if (e.key.code == Keyboard::Key::Up)
				{
					temperature *= 1.0f / 0.9f;
					std::cout << temperature << std::endl;
				}
				else if (e.key.code == Keyboard::Key::Down)
				{
					temperature *= 0.9f;
					std::cout << temperature << std::endl;
				}

				if (e.key.code == Keyboard::Key::Left)
				{
				}
				else if (e.key.code == Keyboard::Key::Right)
				{
					draw_vortices = !draw_vortices;
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
		ImGui::Checkbox("Vortices", &draw_vortices);
		ImGui::End();

		window.clear();

		if (draw_vortices)
			canvas.Draw(squareLattice, squareLattice.getVortices());
		else
			canvas.Draw(squareLattice);

		ImGui::SFML::Render(window);
		window.display();

		machine.Sweep(5000, (double)temperature);

	}

	return 0;
}