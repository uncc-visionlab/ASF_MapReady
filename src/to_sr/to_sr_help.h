#ifndef _TO_SR_HELP_H_
#define _TO_SR_HELP_H_

// NOTES:
//   1. Don't put a '\n' at the end of the string definitions.  Formatting
//      is handled by usage() and print_help()
//   2. To match the formatting found in the other ASF include files and tool
//      output, use the examples below, i.e. a) paragraphs are 3 spaces from
//      the left, b) second line of options starts 13 spaces from the left,
//      c) option descriptions start 8 spaces from the left, etc.
//   3. License, version, and contact info are defined elsewhere (see asf_version.h,
//      svn_rev.h, asf_contact.h, and asf_license.h.  These files are included by
//      help.c
//

// TOOL_NAME is required
#ifdef  TOOL_NAME
#undef  TOOL_NAME
#endif
#define TOOL_NAME \
        "to_sr"

// TOOL_USAGE is required
#ifdef  TOOL_USAGE
#undef  TOOL_USAGE
#endif
#define TOOL_USAGE \
        TOOL_NAME" [-p <pixsize>] <infile> <outfile> [-license] [-version] [-help]"

// TOOL_DESCRIPTION is required
#ifdef  TOOL_DESCRIPTION
#undef  TOOL_DESCRIPTION
#endif
#define TOOL_DESCRIPTION \
    "   Converts an image to slant range.  The input file could be a ground range \n"\
    "   image, or a projected image with state vectors.  Optionally resample the\n"\
    "   output image to a new pixel size.\n"

// TOOL_INPUT is required but is allowed to be an empty string
#ifdef  TOOL_INPUT
#undef  TOOL_INPUT
#endif
#define TOOL_INPUT \
        "   <infile>  The basename of the input data file.  Must be in ASF Internal\n" \
        "        format (.img)."

// TOOL_OUTPUT is required but is allowed to be an empty string
#ifdef  TOOL_OUTPUT
#undef  TOOL_OUTPUT
#endif
#define TOOL_OUTPUT \
        "   <outfile>  The basename of the output file.  The output file will be stored\n" \
        "        in ASF Internal format (.img and .meta)"

// TOOL_OPTIONS is required but is allowed to be an empty string
#ifdef  TOOL_OPTIONS
#undef  TOOL_OPTIONS
#endif
#define TOOL_OPTIONS \
    "   -p <pixel size>\n" \
    "        Specifies the desired slant range pixel size.  If this is not\n"\
    "        specified, it is calculated using the formula:\n"\
    "          (speed of light) / (sample rate * 2*10^6)\n"\
    "   -license\n" \
    "        Print copyright and license for this software then exit.\n" \
    "   -version\n" \
    "        Print version and copyright then exit.\n" \
    "   -help\n" \
    "        Print this help information and then exit."

// TOOL_EXAMPLES is required but is allowed to be an empty string
#ifdef  TOOL_EXAMPLES
#undef  TOOL_EXAMPLES
#endif
#define TOOL_EXAMPLES \
    ""

// TOOL_LIMITATIONS is required but is allowed to be an empty string
#ifdef  TOOL_LIMITATIONS
#undef  TOOL_LIMITATIONS
#endif
#define TOOL_LIMITATIONS \
    ""

// TOOL_SEE_ALSO is required but is allowed to be an empty string
#ifdef  TOOL_SEE_ALSO
#undef  TOOL_SEE_ALSO
#endif
#define TOOL_SEE_ALSO \
    "   sr2gr"

// Prototypes
void check_for_help(int argc, char* argv[]);
void usage();
void print_help();

#endif // _TO_SR_HELP_H_

