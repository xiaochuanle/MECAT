#include "reads_correction_can.h"
#include "reads_correction_m4.h"

int main(int argc, char** argv)
{
    ReadsCorrectionOptions rco;
	parse_reads_correction_arguments(argc, argv, rco);
	if (rco.input_type == INPUT_TYPE_CAN)
	{
		return reads_correction_can(rco);
	}
	else
	{
		return reads_correction_m4(rco);
	}
	return 0;
}
