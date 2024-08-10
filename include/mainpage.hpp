// TODO: finish this page
/** @mainpage
 *
 * @authors Bin-Hui Chen (orcid: https://orcid.org/0000-0001-8962-663X)
 *
 * @section intro Introduction
 * This package provides code templates for use by GLAST developers.
 * All header files for external access are located in the templates directory,
 * as it is customary to put external public header files in the packageName
 * directory. Header files that are not meant for external access reside in
 * the src directory. Source files are located in the src directory. Files
 * related to loading sharable libraries are located in the src/Dll directory.
 * There are 3 examples:
 *
 * - User-Defined generic C++ class
 * -# templates/ClassTemplate.h
 * -# src/ClassTemplate.cxx
 * - User-Defined Gaudi Algorithm
 * -# src/ExampleAlg.cxx
 * -# src/Dll/templates_dll.cxx
 * -# src/Dll/templates_load.cxx
 * - User-Defined Gaudi Service
 * -# templates/IExampleSvc.h
 * -# templates/ExampleSvc.h
 * -# src/ExampleSvc.cxx
 * -# src/Dll/templates_dll.cxx
 * -# src/Dll/templates_load.cxx
 *
 *
 * Also note the existence of the following directories:
 * - cmt
 * -# Contains the requirements file
 * - doc
 * -# Contains the release.notes file
 *
 *
 * As you prepare to develop code for GLAST SAS, please be sure you are aware
 * of our current
 * <A HREF="http://www-glast.slac.stanford.edu/software/CodeHowTo/codeStandards.html"> Coding
 * Standards </A>
 *
 *
 * If using the code in this package as an example - please modify the comments
 * as appropriate for your own specific code.
 *
 * <hr>
 * @section notes release.notes
 * release.notes
 * <hr>
 * @section requirements requirements
 * @verbinclude requirements
 * <hr>
 * @todo [optionally include text about more work to be done]
 * @todo Give each todo item its own line
 *
 */
