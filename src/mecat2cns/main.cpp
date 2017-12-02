#include "reads_correction_can.h"
#include "reads_correction_m4.h"

int main(int argc, char** argv)
{
    ReadsCorrectionOptions rco;
	int r = parse_arguments(argc, argv, rco);
	print_options(rco);
	if (r) {
		print_usage(argv[0]);
		exit(1);
	}
	if (rco.print_usage_info) {
		print_usage(argv[0]);
		exit(0);
	}
	
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
