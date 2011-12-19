/*=========================================================================

Program:   EpilepsyViewerInteraction
Module:    EpilepsyViewerInteraction.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

#include "EpilepsyViewerInteraction.h"

// from ToolCursor
#include "vtkToolCursor.h"
#include "vtkWindowLevelTool.h"
#include "vtkSliceImageTool.h"
#include "vtkRotateCameraTool.h"
#include "vtkPanCameraTool.h"
#include "vtkZoomCameraTool.h"

//----------------------------------------------------------------------------
EpilepsyViewerInteraction::EpilepsyViewerInteraction()
{
  this->ToolCursor = vtkSmartPointer<vtkToolCursor>::New();
  this->ToolCursor->SetScale(1.0);
}

//----------------------------------------------------------------------------
EpilepsyViewerInteraction::~EpilepsyViewerInteraction()
{
}

//----------------------------------------------------------------------------
vtkToolCursor *EpilepsyViewerInteraction::GetToolCursor()
{
  return this->ToolCursor;
}

//----------------------------------------------------------------------------
void EpilepsyViewerInteraction::BindInteraction(vtkRenderer *renderer)
{
  vtkToolCursor *cursor = this->ToolCursor;

  cursor->SetRenderer(renderer);

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
}
