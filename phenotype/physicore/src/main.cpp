#include <microenvironment.h>

using namespace physicore;

int main(int argc, char** argv)
{
	(void)argv;
	(void)argc;

	biofvm::microenvironment m(0.1);

	m.run_single_timestep();

	return 0;
}
