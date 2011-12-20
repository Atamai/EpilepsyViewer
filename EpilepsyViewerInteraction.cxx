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
#include "vtkPushPlaneTool.h"
#include "vtkRotateCameraTool.h"
#include "vtkPanCameraTool.h"
#include "vtkZoomCameraTool.h"
#include "vtkActionCursorShapes.h"
#include "vtkGeometricCursorShapes.h"

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

  // cursor shapes
  vtkSmartPointer<vtkActionCursorShapes> actionShapes =
    vtkSmartPointer<vtkActionCursorShapes>::New();
  int rotateShape = cursor->AddShape(actionShapes, "Rotate");
  int pusherShape = cursor->AddShape(actionShapes, "Push");
  int zoomShape = cursor->AddShape(actionShapes, "Zoom");
  int panShape = cursor->AddShape(actionShapes, "Move");

  vtkSmartPointer<vtkGeometricCursorShapes> geometricShapes =
    vtkSmartPointer<vtkGeometricCursorShapes>::New();
  int coneShape = cursor->AddShape(geometricShapes, "Cone");
  int crossShape = cursor->AddShape(geometricShapes, "SplitCross");

  // tools
  vtkSmartPointer<vtkPushPlaneTool> pushTool =
    vtkSmartPointer<vtkPushPlaneTool>::New();
  int pushId = cursor->AddAction(pushTool);

  vtkSmartPointer<vtkRotateCameraTool> rotateTool =
    vtkSmartPointer<vtkRotateCameraTool>::New();
  int rotateId = cursor->AddAction(rotateTool);

  vtkSmartPointer<vtkPanCameraTool> panTool =
    vtkSmartPointer<vtkPanCameraTool>::New();
  int panId = cursor->AddAction(panTool);

  vtkSmartPointer<vtkZoomCameraTool> zoomTool =
    vtkSmartPointer<vtkZoomCameraTool>::New();
  int zoomId = cursor->AddAction(zoomTool);

  //vtkSmartPointer<vtkWindowLevelTool> winlevTool =
  //  vtkSmartPointer<vtkWindowLevelTool>::New();
  //int winlevId = cursor->AddAction(winlevTool);

  // context
  int actorInfo = 0; // what is under the cursor?
  int modifier = 0; // what modifier keys?
  int mode = 0; // what is the current mode?

  // ------------ Camera Bindings-------------
  actorInfo = 0;

  // rotation
  modifier = 0;
  cursor->BindShape(rotateShape, mode, actorInfo, modifier | VTK_TOOL_B1);
  cursor->BindAction(rotateId, mode, actorInfo, modifier | VTK_TOOL_B1);

  // panning
  modifier = VTK_TOOL_SHIFT;
  cursor->BindShape(panShape, mode, actorInfo, modifier);
  cursor->BindShape(panShape, mode, actorInfo, modifier | VTK_TOOL_B1);
  cursor->BindAction(panId, mode, actorInfo, modifier | VTK_TOOL_B1);

  // Also bind to middle button
  modifier = 0;
  cursor->BindShape(panShape, mode, actorInfo, modifier | VTK_TOOL_B3);
  cursor->BindAction(panId, mode, actorInfo, modifier | VTK_TOOL_B3);

  // Binding for "Zoom" cursor and action
  modifier = VTK_TOOL_CONTROL;
  cursor->BindShape(zoomShape, mode, actorInfo, modifier);
  cursor->BindShape(zoomShape, mode, actorInfo, modifier | VTK_TOOL_B1);
  cursor->BindAction(zoomId, mode, actorInfo, modifier | VTK_TOOL_B1);

  // Also bind to right button
  modifier = 0;
  cursor->BindShape(zoomShape, mode, actorInfo, modifier | VTK_TOOL_B2);
  cursor->BindAction(zoomId, mode, actorInfo, modifier | VTK_TOOL_B2);

  // ------------ Default Prop3D Bindings-------------
  actorInfo = VTK_TOOL_PROP3D;
  modifier = 0;
  cursor->BindShape(coneShape, mode, actorInfo, modifier);
  cursor->BindShape(rotateShape, mode, actorInfo, modifier | VTK_TOOL_B1);
  cursor->BindAction(rotateId, mode, actorInfo, modifier | VTK_TOOL_B1);

  // ------------ Bindings for Volume -------------
  actorInfo = (VTK_TOOL_VOLUME | VTK_TOOL_CLIP_PLANE | VTK_TOOL_CROP_PLANE);
  modifier = 0;
  cursor->BindShape(pusherShape, mode, actorInfo, modifier);
  cursor->BindShape(pusherShape, mode, actorInfo, modifier | VTK_TOOL_B1);
  cursor->BindAction(pushId, mode, actorInfo, modifier | VTK_TOOL_B1);

  // ------------ Bindings for ImageActor -------------
  actorInfo = VTK_TOOL_IMAGE_ACTOR;
  modifier = 0;
  cursor->BindShape(pusherShape, mode, actorInfo, modifier);
  cursor->BindShape(pusherShape, mode, actorInfo, modifier | VTK_TOOL_B1);
  cursor->BindAction(pushId, mode, actorInfo, modifier | VTK_TOOL_B1);
}
