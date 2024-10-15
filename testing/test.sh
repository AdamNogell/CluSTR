#!/bin/bash

function usage()
{
	echo "CluSTR Tool"
	echo "Analysis of Indels in Short Tandem Repeats"
	echo
	echo "Usage: $0 [-i INPUT.TXT] [-t TITLE] [-p] [-o OUTPUT]"
	echo
	echo "  [-h]           - Display this help message"
    echo "  [-v]           - Show version information"
	echo "  [-i INPUT.TXT] - Input file. Should be the output of \"solexa_quality_statistics\" program."
	echo "  [-o OUTPUT]    - Output file name. default is STDOUT."
	echo "  [-t TITLE]     - Title (usually the solexa file name) - will be plotted on the graph."
	echo
	exit 
}
if [ $# -eq 0 ]; then
    usage
    exit 1
fi

while getopts ":hvf:" opt; do
    case ${opt} in
        h )
            usage
            exit 0
            ;;
        v )
            echo "CluSTR Tool"
            echo "Version 0.0"
            exit 0
            ;;
        f )
            echo "Processing file: $OPTARG"
            ;;
        \? )
            printf "\e[31mInvalid option:\e[0m -%s\n" "$OPTARG" >&2
            echo
            echo "Usage: $0 [-i INPUT.TXT] [-t TITLE] [-p] [-o OUTPUT]"
            echo
            echo "Use [-h] to display help message"
            echo -e "\e[31mThis is red text\e[0m"
            exit 1
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done