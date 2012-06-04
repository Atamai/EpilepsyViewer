/*=========================================================================

Program:   EpilepsyViewerDisplay
Module:    EpilepsyViewerDisplay.h

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

#ifndef __EpilepsyViewerDisplay_h
#define __EpilepsyViewerDisplay_h

#include <vtkSmartPointer.h>
#include <vector>

class vtkImageReslice;
class vtkLODProp3D;
class vtkVolumeProperty;
class vtkActor;
class vtkProperty;
class vtkPlane;
class vtkImageData;
class vtkImageStack;
class vtkImageSlice;
class vtkImageProperty;
class vtkImageResliceMapper;
class vtkRenderWindow;
class vtkRenderer;
class vtkPlaneCollection;
class vtkMatrix4x4;

class EpilepsyViewerData;

//! Handle all data display for the EpilepsyViewer
class EpilepsyViewerDisplay
{
public:
  EpilepsyViewerDisplay();
  ~EpilepsyViewerDisplay();

  //! The three slice orientations.
  enum Orientation { Sagittal, Coronal, Axial, OrientationCount };

  //! The layers.
  enum Layer { MRLayer, CTLayer };

  //! Set the render window.  This can only be called once.
  void SetRenderWindow(vtkRenderWindow *renwin);

  //! Get the the main renderer.
  vtkRenderer *GetMainRenderer();

  //! Set the data that is to be viewed.
  void SetData(EpilepsyViewerData *data);

  //! Set the orientation of all slices from the provided matrix.
  void SetSliceOrientation(vtkMatrix4x4 *matrix);

  //! Set the orientation of just a single slice.
  void SetSliceOrientation(int orientation, vtkMatrix4x4 *matrix);

  //! This method is called on every render.
  void OnRender();

  //! Get the CT image display properties object
  vtkImageProperty *GetCTImageProperty();

  //! Get the MR image display properties object
  vtkImageProperty *GetMRImageProperty();

  //! Get the MR volume display properties object
  vtkVolumeProperty *GetMRVolumeProperty();

  //! Get the electrode display properties object
  vtkProperty *GetElectrodesProperty();

private:
  //! Add a layer to the viewer, using the provided property
  //! which should have a LayerNumber set.
  void AddLayer(vtkImageProperty *property);

  //! Get the vtkImageSlice object for the specified layer and orientation.
  vtkImageSlice *GetImageSlice(int layer, int orientation);

  //! Generate a set of clipping planes from a set of bounds and a
  //! transformation matrix.
  static void GenerateClippingPlanes(vtkPlaneCollection *planes,
    vtkMatrix4x4 *matrix, const double bounds[6]);

  //! Tighten the clipping planes around a volume so that they touch on
  //! each side.
  static void TightenClippingPlanes(
    vtkPlaneCollection *planes, vtkMatrix4x4 *matrix, vtkImageData *data);

  //! Get a transformed set of bounds that encloses the original bounds.
  static void TransformBounds(vtkMatrix4x4 *matrix,
    const double inbounds[6], double outbounds[6]);

  //! Set up vtkImageReslice so that it extracts the specified plane.
  static void SetReslicePlane(
    vtkImageReslice *reslice, const double plane[4]);

  //! Convert a plane to a reslice matrix.
  static void ConvertPlaneToResliceAxes(
    const double plane[4], double matrix[16]);

  vtkSmartPointer<vtkRenderWindow> RenderWindow;
  vtkSmartPointer<vtkRenderer> MainRenderer;

  vtkSmartPointer<vtkImageProperty> MRImageProperty;
  vtkSmartPointer<vtkImageProperty> CTImageProperty;
  vtkSmartPointer<vtkVolumeProperty> MRVolumeProperty;
  vtkSmartPointer<vtkProperty> ElectrodesProperty;

  vtkSmartPointer<vtkImageReslice> BrainVolumeReslice;
  vtkSmartPointer<vtkLODProp3D> BrainVolume;
  std::vector<int> BrainVolumeLODIds;
  vtkSmartPointer<vtkImageReslice> HeadVolumeReslice;
  vtkSmartPointer<vtkLODProp3D> HeadVolume;
  std::vector<int> HeadVolumeLODIds;

  vtkSmartPointer<vtkActor> ElectrodesActor;

  //! A struct to hold everything for an ortho plane slice.
  struct Slice
  {
    Slice(int orientation);
    virtual ~Slice();

    vtkSmartPointer<vtkImageStack> Stack;
    vtkSmartPointer<vtkPlane> Plane;
    int Orientation;
  };

  //! A vector containing all of the slices.
  std::vector<Slice> Slices;
};

#endif
