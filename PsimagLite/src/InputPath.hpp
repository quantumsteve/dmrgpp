#ifndef INPUTPATH_HPP
#define INPUTPATH_HPP
#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

namespace PsimagLite {

class InputPath {
public:

	InputPath()
	    : paths_(1, "")
	{ }

	void push(const std::string& path)
	{
		// path already in paths
		if (std::find(paths_.begin(), paths_.end(), path) != paths_.end()) {
			return;
		}

		paths_.push_back(path);
	}

	std::string findFirst(const std::string& file) const
	{
		for (auto rit = paths_.rbegin(); rit != paths_.rend(); ++rit) {
			const std::string&    dir       = *rit;
			std::string           test_file = (dir.empty()) ? file : dir + "/" + file;
			std::filesystem::path path      = test_file;
			if (std::filesystem::exists(path)) {
				return test_file;
			}
		}

		throw std::runtime_error("File " + file
		                         + " could not be found within available paths\n");
	}

private:

	std::vector<std::string> paths_;
};
}
#endif // INPUTPATH_HPP
