// ptx2png.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "lodepng.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#define M_PI 3.1415926535


int main(int argc, char* argv[])
{

	if (argc < 2) {
		cout << "Usage: ptx2png infile outfile\n  infile: the input file (*.ptx)\n  outfile: the output file (*.png)" << endl;
		exit(1);
	}

	ifstream ifs;
	ifs.open(argv[1]);

	if (ifs.is_open()) {
		cout << "Converting " << argv[1] << " to image " << argv[2] << endl;

		string line;
		getline(ifs, line);
		const int width = stoi(line);
		getline(ifs, line);
		const int height = stoi(line);

		cout << "Image dimensions: " << width <<  " x " << height << endl;
		
		// skips header lines
		for(int i=0; i<8; i++) {
			getline(ifs, line);
		}

		unsigned int total_lines = width * height;
		vector<unsigned char> buf(4*total_lines);

		double cx, cy, cz, intensity;
		int r, g, b;
		for (int x = 0; x < width; ++x) {
			cout << '\r' << x << '/'<< width;
			for (int y = height-1; y >= 0; --y) {

				ifs >> cx >> cy >> cz >> intensity >> r >> g >> b;

				// need to compute index because ptx is ordered by columns and png is ordered by rows
				int index = 4*(y*width+x);

				buf[index] = static_cast<unsigned char>(r);
				buf[index+1] = static_cast<unsigned char>(g);
				buf[index+2] = static_cast<unsigned char>(b);
				buf[index+3] = 255;
			}
		}
		cout << endl;

		vector<unsigned char> png;
		lodepng::State state;
		unsigned int error = lodepng::encode(png, buf, width, height, state);

		if (!error) {
			lodepng::save_file(png, argv[2]);
			cout << "SUCCESS: saved image file " << argv[2] << endl;
		} else {
			cout << "ERROR: error while encoding, this may be due to an incorrect specified number of pixels" << endl;
		}

		ifs.close();
	}

	return 0;
}



