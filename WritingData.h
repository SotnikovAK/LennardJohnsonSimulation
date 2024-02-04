#pragma once

#include "Libraries.h"
#include "Element.h"


class WritingData {
private:
	std::ofstream outputFile;

public:
	WritingData(const std::string& filename) {
		outputFile.open(filename, std::ios_base::app);
	}
	void printing(const std::string& message) {
		if (outputFile.is_open()) {
			outputFile << message << std::endl;

		}

	}

	~WritingData() {
		if (outputFile.is_open()) {
			outputFile.close();

		}
	}
};

class SpecialWritingFunction {

private:

public:
	SpecialWritingFunction() {};

	void PrintEnergy(int t, float U, std::list<ElementarElement> ElementList);

	void PrintDiffusion(int t, std::list<ElementarElement> ElementList);

	void PrintV(int t, std::list<ElementarElement> ElementList);

	void PrintImp(int t, std::list<ElementarElement> ElementList);


};

