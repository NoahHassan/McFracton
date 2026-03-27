#include <SFML/Graphics.hpp>
#include <iostream>

#include "Canvas.h"
#include "XYSquare.h"
#include "Timer.h"
#include "McMachine.h"

int main() {

	using namespace sf;

	RenderWindow window(VideoMode(1200, 800), "Simulation");
	window.setVerticalSyncEnabled(true);

	XYSquare squareLattice(80);
	Canvas canvas(window);
	canvas.Initialize(squareLattice, 7.0);

	McMachine::NumericalParams params(
		100.0, 0.1, 0.1
	);
	McMachine machine(params, squareLattice);

	bool draw_vortices = true;

	Timer timer;
	double temperature = 100.0;
	while (window.isOpen())
	{
		Event e;
		while (window.pollEvent(e))
		{
			if (e.type == Event::KeyPressed)
			{
				if (e.key.code == Keyboard::Key::Enter)
				{
					std::cout << squareLattice.getEnergy() << std::endl;
				}

				if (e.key.code == Keyboard::Key::Up)
				{
					temperature *= 1.0 / 0.9;
					std::cout << temperature << std::endl;
				}
				else if (e.key.code == Keyboard::Key::Down)
				{
					temperature *= 0.9;
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

		machine.Sweep(5000, temperature);

		if (timer.elapsed() >= 1.0f / 30.0f)
		{
			if (draw_vortices)
				canvas.Draw(squareLattice, squareLattice.getVortices());
			else
				canvas.Draw(squareLattice);
			timer.reset();
		}

	}

	return 0;

	/*
	Pseudo-Code for how the Code should work:

	define Lattice
	Lattice contains:
	  - site fields (exposed as index)
	  - plaquette fields (exposed as index)
	  - float proposeFlip() (returns energy difference)
	  - getEnergy() (returns complete energy)
	  - internally has a function or a struct mapping a site/plaquette
	    to connected objects that will be affected

	initialize McMachine(Params, Lattice)
	McMachine contains:
	  - Params struct (t_max, t_min, sweeps, etc.)
	  - internally: proposes Flips and updates them

	optional drawing stuff
	draw(Lattice)
	*/
}