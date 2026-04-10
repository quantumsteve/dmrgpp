#include "InputPath.hpp"
#include <iostream>
#include <unistd.h>

int main(int argc, char* argv[])
{
	PsimagLite::InputPath input_path;
	std::string           usage
	    = std::string("USAGE: ") + argv[0] + " -I path1 -I path2 ... -f filename\n";

	std::string filename;
	int         opt = 0;
	while ((opt = getopt(argc, argv, "f:I:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'I':
			input_path.push(optarg);
			break;
		default:
			throw std::runtime_error(usage);
		}
	}

	if (filename.empty()) {
		throw std::runtime_error("File must be provided\n");
	}

	std::string result = input_path.findFirst(filename);
	std::cout << result << "\n";
}
