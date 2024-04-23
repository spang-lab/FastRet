exitcode=0
red="\033[0;31m"
green="\033[0;32m"
reset="\033[0m"
ok="${green}ok${reset}"
failed="${red}failed${reset}"

if [ -z "$1" ]
then
    echo -e "${red}Error: no version provided${reset}"
    echo -e "${red}Usage: install_jdk (default-jdk|openjdk-8-jdk) [--verbose]${reset}"
    echo -e "${red}Valid versions: 8, 11, 17, 18, 19, 21${reset}"
    exit 1
fi

valid_versions=("default-jdk" "openjdk-8-jdk" "openjdk-11-jdk" "openjdk-17-jdk" "openjdk-18-jdk" "openjdk-19-jdk" "openjdk-21-jdk")
if [[ ! " ${valid_versions[@]} " =~ " $1 " ]]; then
    echo -e "${red}Error: Invalid JDK version${reset}"
    echo -e "${red}Valid versions are: ${valid_versions[@]}${reset}"
    exit 1
fi

if [[ "$@" == *"--verbose"* ]]; then
    echo "Verbose mode enabled"
    OUTLOG="/dev/stdout"
    ERRLOG="/dev/stdout"
    LINEEND="\n"
else
    OUTLOG="/dev/null"
    ERRLOG="/dev/null"
    LINEEND=""
fi

echo -n "Removing old java installations ... $LINEEND"
apt-get update
apt-get remove -y '*j[rd]k*' 1>$OUTLOG 2>$ERRLOG
command -v java 1>$OUTLOG 2>$ERRLOG
if [ $? -ne 0 ]; then echo -e $ok; else echo -e $failed; exitcode=1; fi

echo -n "Installing $1 ... $LINEEND"
apt-get install --no-install-recommends -y "$1" 1>$OUTLOG 2>$ERRLOG
command -v java 1>$OUTLOG 2>$ERRLOG
if [ $? -eq 0 ]; then echo -e $ok; else echo -e $failed; exitcode=1; fi

echo -n "Calling 'R CMD javareconf -e' ... $LINEEND"
R CMD javareconf -e 1>$OUTLOG 2>$ERRLOG
if [ $? -eq 0 ]; then echo -e $ok; else echo -e $failed; exitcode=1; fi

exit $exitcode
