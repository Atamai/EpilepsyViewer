/*=========================================================================

Program:   EpilepsyViewerInteraction
Module:    EpilepsyViewerInteraction.h

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

#ifndef __EpilepsyViewerInteraction_h
#define __EpilepsyViewerInteraction_h

// from VTK
#include <vtkSmartPointer.h>

class vtkToolCursor;
class vtkRenderer;

//! Handle all user interaction with the render window
class EpilepsyViewerInteraction
{
public:
  EpilepsyViewerInteraction();
  ~EpilepsyViewerInteraction();

  //! Get the tool cursor that controls the interaction.
  vtkToolCursor *GetToolCursor();

  //! Bind the interaction.
  void BindInteraction(vtkRenderer *renderer);

private:
  vtkSmartPointer<vtkToolCursor> ToolCursor;
};

#endif
