#include "WritingData.h"

void SpecialWritingFunction::PrintEnergy(int t, float U, std::list<ElementarElement> ElementList)
{
	WritingData file("./DATA/Energy.txt");

	float SumKinEnergy = 0.0;

	for (auto element = ElementList.begin(); element != ElementList.end(); element++)
	{
		SumKinEnergy += element->KinEnergy();
	}
	file.printing(std::to_string(t) + " " + std::to_string(SumKinEnergy) + " " + std::to_string(U));

	WritingData file2("./DATA/FullEnergy.txt");
	file2.printing(std::to_string(t) + " " + std::to_string(SumKinEnergy + U));
}

void SpecialWritingFunction::PrintDiffusion(int t, std::list<ElementarElement> ElementList)
{
	WritingData file("./DATA/Dispersia.txt");

	float SumDiff = 0.0;
	for (auto element = ElementList.begin(); element != ElementList.end(); element++)
	{
		SumDiff += element->Disp();
	}
	SumDiff /= ElementList.size();
	file.printing(std::to_string(t) + " " + std::to_string(SumDiff) );

}

void SpecialWritingFunction::PrintV(int t, std::list<ElementarElement> ElementList)
{
	WritingData file("./DATA/AbsV.txt");

	for (auto element = ElementList.begin(); element != ElementList.end(); element++)
	{
		file.printing(std::to_string(element->absV()));
		
	}
}

void SpecialWritingFunction::PrintImp(int t, std::list<ElementarElement> ElementList)
{
	WritingData file("./DATA/P.txt");

	float SumP[3] = {0.0,0.0,0.0};
	for (auto element = ElementList.begin(); element != ElementList.end(); element++)
	{
		float P[3] = { 0.0,0.0,0.0 };

		element->Imp(P);

		std::tie(SumP[0], SumP[1], SumP[2]) = std::make_tuple(SumP[0] + P[0], SumP[1] + P[1], SumP[2] + P[2]);
	}
	file.printing(std::to_string(t) +  " " + std::to_string(SumP[0]) + " " + std::to_string(SumP[1]) + " " + std::to_string(SumP[2]));
}



