/*=========================================================================

Program:   EpilepsyViewer
Module:    EpilepsyViewer.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

#include "EpilepsyViewerData.h"
#include "EpilepsyViewerDisplay.h"
#include "EpilepsyViewerInteraction.h"

// from VTK
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

// from ToolCursor
#include "vtkToolCursorInteractorObserver.h"

int main (int argc, char *argv[])
{
  if (argc < 2)
    {
    cout << "Usage: " << argv[0] << " <EpilepsyDirectory>" << endl;
    return EXIT_FAILURE;
    }

  // data, display and interaction managers
  EpilepsyViewerData dataManagerObject;
  EpilepsyViewerDisplay displayManagerObject;
  EpilepsyViewerInteraction interactionManagerObject;

  // get pointers to the objects (for stylistic reasons only)
  EpilepsyViewerData *dataManager = &dataManagerObject;
  EpilepsyViewerDisplay *displayManager = &displayManagerObject;
  EpilepsyViewerInteraction *interactionManager = &interactionManagerObject;

  // set up the data first
  dataManager->SetDataDirectory(argv[1]);

  if (!dataManager->LoadFromDataDirectory())
    {
    cerr << "Unable to read data, exiting.\n";
    return EXIT_FAILURE;
    }

  // go through the data processing steps
  dataManager->RegisterMRHeadToCTHead();
  dataManager->ExtractMRBrain();
  dataManager->ExtractCTElectrodes();

  // push the data to the display manager
  displayManager->SetData(dataManager);

  // create the render window
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> interactor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindow->SetInteractor(interactor);

  // for "basic" epilepsy viewer, use interactor observer
  vtkSmartPointer<vtkToolCursorInteractorObserver> observer =
    vtkSmartPointer<vtkToolCursorInteractorObserver>::New();
  observer->SetToolCursor(interactionManager->GetToolCursor());
  observer->SetInteractor(interactor);
  observer->SetEnabled(1);

  // put it all together and display it
  renderWindow->SetSize(1024,1024);
  displayManager->SetRenderWindow(renderWindow);
  interactionManager->BindInteraction(displayManager->GetMainRenderer());
  renderWindow->Render();
  interactor->Start();

  return 0;
}
