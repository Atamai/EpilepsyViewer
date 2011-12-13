/*=========================================================================

Program:   EpilepsyViewer
Module:    EpilepsyViewer.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

#include "EpilepsyViewerData.h"
#include "EpilepsyViewerDisplay.h"

// from ToolCursor
#include "vtkToolCursor.h"
#include "vtkWindowLevelTool.h"
#include "vtkSliceImageTool.h"
#include "vtkRotateCameraTool.h"
#include "vtkPanCameraTool.h"
#include "vtkZoomCameraTool.h"
#include "vtkToolCursorInteractorObserver.h"

// from VTK
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

int main (int argc, char *argv[])
{
  if (argc < 2)
    {
    cout << "Usage: " << argv[0] << " <EpilepsyDirectory>" << endl;
    return EXIT_FAILURE;
    }

  // data manager and display manager
  EpilepsyViewerData dataManagerObject;
  EpilepsyViewerDisplay displayManagerObject;
  EpilepsyViewerData *dataManager = &dataManagerObject;
  EpilepsyViewerDisplay *displayManager = &displayManagerObject;

  dataManager->SetDataDirectory(argv[1]);

  if (!dataManager->LoadFromDataDirectory())
    {
    cerr << "Unable to read data, exiting.\n";
    return EXIT_FAILURE;
    }

  dataManager->RegisterMRHeadToCTHead();
  dataManager->ExtractMRBrain();
  dataManager->ExtractCTElectrodes();
  displayManager->SetData(dataManager);

  // add ToolCursor items
  vtkSmartPointer<vtkToolCursor> cursor =
    vtkSmartPointer<vtkToolCursor>::New();
  cursor->SetRenderer(displayManager->GetMainRenderer());
  cursor->SetScale(1.0);

  vtkSmartPointer<vtkWindowLevelTool> winlevTool =
    vtkSmartPointer<vtkWindowLevelTool>::New();
  int winlevId = cursor->AddAction(winlevTool);

  vtkSmartPointer<vtkSliceImageTool> sliceTool =
    vtkSmartPointer<vtkSliceImageTool>::New();
  sliceTool->JumpToNearestSliceOn();
  int sliceId = cursor->AddAction(sliceTool);

  vtkSmartPointer<vtkRotateCameraTool> rotateTool =
    vtkSmartPointer<vtkRotateCameraTool>::New();
  int rotateId = cursor->AddAction(rotateTool);

  vtkSmartPointer<vtkPanCameraTool> panTool =
    vtkSmartPointer<vtkPanCameraTool>::New();
  int panId = cursor->AddAction(panTool);

  vtkSmartPointer<vtkZoomCameraTool> zoomTool =
    vtkSmartPointer<vtkZoomCameraTool>::New();
  int zoomId = cursor->AddAction(zoomTool);

  // Bind all the tools
  cursor->BindDefaultActions();
/*
  cursor->BindAction(rotateId, 0, 0, VTK_TOOL_B1);
  cursor->BindAction(sliceId, 0, 0, VTK_TOOL_SHIFT | VTK_TOOL_B1);
  cursor->BindAction(panId, 0, 0, VTK_TOOL_SHIFT | VTK_TOOL_B2);
  cursor->BindAction(zoomId, 0, 0, VTK_TOOL_B2);
*/

  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> interactor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindow->SetInteractor(interactor);

  renderWindow->SetSize(1024,1024);

  vtkSmartPointer<vtkToolCursorInteractorObserver> observer =
    vtkSmartPointer<vtkToolCursorInteractorObserver>::New();
  observer->SetToolCursor(cursor);
  observer->SetInteractor(interactor);
  observer->SetEnabled(1);

  displayManager->SetRenderWindow(renderWindow);
  renderWindow->Render();
  interactor->Start();

  return 0;
}
