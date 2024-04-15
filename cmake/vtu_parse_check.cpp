/**
 * Check that VTK can successfully parse a .vtu file containing binary data, in case broken by a libexpat1 patch.
 * VTK prints XML parsing errors to stdout by default, with no programmatically accessible way to detect that one has occurred (?!) without setting a custom errorObservers?
*/

// #include <vtkUnstructuredGrid.h>
// #include <vtkUnstructuredGridReader.h>
#include <vtkCommand.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <cstdio>
#include <string>

/**
 * VTK custom error observer adapted from VTK example https://examples.vtk.org/site/Cxx/Utilities/ObserveError/
*/
class ErrorObserver : public vtkCommand
{
 public:
    ErrorObserver()
    : Error(false), Warning(false), ErrorMessage(""), WarningMessage("")
    {
    }
    static ErrorObserver* New() { return new ErrorObserver; }
    bool GetError() const { return this->Error; }
    bool GetWarning() const { return this->Warning; }
    void Clear()
    {
        this->Error = false;
        this->Warning = false;
        this->ErrorMessage = "";
        this->WarningMessage = "";
    }
    virtual void Execute(vtkObject* vtkNotUsed(caller), unsigned long event,
    void* calldata)
    {
        switch (event)
        {
            case vtkCommand::ErrorEvent:
                ErrorMessage = static_cast<char*>(calldata);
                this->Error = true;
                break;
            case vtkCommand::WarningEvent:
                WarningMessage = static_cast<char*>(calldata);
                this->Warning = true;
                break;
        }
    }
    std::string GetErrorMessage() { return ErrorMessage; }
    std::string GetWarningMessage() { return WarningMessage; }

 private:
    bool Error;
    bool Warning;
    std::string ErrorMessage;
    std::string WarningMessage;
};

/**
 * Main method taking a path to a .vtu file as an argument. Return code indicates if VTU parsing resulted in any errors or not.
 */
int main(int argc, char * argv[]) {
    // Ensure an argument is provided
    if (argc != 2) {
        fprintf(stderr, "Error: a single argument is required (path to .vtu file)\n");
        return EXIT_FAILURE;
    }
    // Declare the xml reader
    vtkXMLUnstructuredGridReader* reader = vtkXMLUnstructuredGridReader::New();
    // Declare and attach custom VTK error observers for xml reading and parsing
    ErrorObserver* readerErrorObserver = ErrorObserver::New();
    reader->SetReaderErrorObserver(readerErrorObserver);
    ErrorObserver* parserErrorObserver = ErrorObserver::New();
    reader->SetParserErrorObserver(parserErrorObserver);
    // Set the filename 
    reader->SetFileName(argv[1]);
    // Parse the .vtu file
    reader->Update();
    // Error of any xml reading errors occurred
    if (readerErrorObserver->GetError())
    {
        return EXIT_FAILURE;
    }
    // Error if any xml parsing errors occurred (i.e. patched libexpat1 without patched vtk)
    if (parserErrorObserver->GetError())
    {
        return EXIT_FAILURE;
    }
    // If no errors detected, return success
    return EXIT_SUCCESS;
}