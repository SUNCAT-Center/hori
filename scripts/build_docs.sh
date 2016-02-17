####################################################################################################
# RUN SPHINX APIDOC AND BUILD
# - requires that DEV_PATH is an environment variable
#
# USAGE:
# - build_docs.sh [options]
#
# OPTIONS:
# -o
#	opens index.html and output.txt immediately after building
####################################################################################################


#...................................................................................................
# CONFIGURATION
#...................................................................................................
SOURCE_DIR="hori/"
DOC_PATH="${DEV_PATH}/${SOURCE_DIR}/docs"
SOURCE_PATH="${DOC_PATH}/source"
HTML_PATH="${DOC_PATH}/html"


#...................................................................................................
# SETUP
#...................................................................................................
mkdir -p $DOC_PATH
cd $DOC_PATH
echo ">>> Cleaning build directory..."
rm -rf ${SOURCE_PATH}/_dynamic ${HTML_PATH}/*


#...................................................................................................
# BUILD
#...................................................................................................
echo ">>> Writing source documentation files..."
sphinx-apidoc -M -o ${SOURCE_PATH}/_dynamic ${DEV_PATH}/${SOURCE_DIR} > /dev/null # docs in _dynamic

echo ">>> Running doctests..."
sphinx-build -b doctest -Q ${SOURCE_PATH} /tmp/doctest
tail -6 /tmp/doctest/output.txt
echo ">>> Doctest run complete. Full doctest output in /tmp/doctest/output.txt..."

echo ">>> Building html documentation files..."
sphinx-build -b html -Q ${SOURCE_PATH} ${HTML_PATH} # cleaner console output than make html
echo ">>> Build complete. The HTML pages are in ${HTML_PATH}..."


#...................................................................................................
# WRAPUP
#...................................................................................................
while getopts ":o" opt; do
	case $opt in
		o) open ${HTML_PATH}/index.html /tmp/doctest/output.txt ;;
	esac
done
git add ${HTML_PATH}
cd - > /dev/null
