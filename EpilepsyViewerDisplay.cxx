/*=========================================================================

Program:   EpilepsyViewerDisplay
Module:    EpilepsyViewerDisplay.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

#include "EpilepsyViewerDisplay.h"
#include "EpilepsyViewerData.h"

// From VTK
#include <vtkImageReslice.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkLODProp3D.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkVolumeTextureMapper3D.h>
#include <vtkVolumeProperty.h>
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkImageStack.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceCollection.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageProperty.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkPlane.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkMath.h>
#include <vtkCommand.h>

// From standard library
#include <stddef.h>

//----------------------------------------------------------------------------
EpilepsyViewerDisplay::Slice::Slice(int orientation)
{
  this->Stack = vtkSmartPointer<vtkImageStack>::New();
  this->Plane = vtkSmartPointer<vtkPlane>::New();
  this->Orientation = orientation;
  if (orientation >= 0 && orientation < 2)
    {
    double normal[3];
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    normal[orientation] = 1.0;
    this->Plane->SetNormal(normal);
    }
}

//----------------------------------------------------------------------------
EpilepsyViewerDisplay::Slice::~Slice()
{
}

//----------------------------------------------------------------------------
EpilepsyViewerDisplay::EpilepsyViewerDisplay()
{
  // the renderers
  this->MainRenderer = vtkSmartPointer<vtkRenderer>::New();
  this->MainRenderer->SetBackground(0.2, 0.2, 0.2);

  // all of the properties that the App has access to
  this->MRImageProperty = vtkSmartPointer<vtkImageProperty>::New();
  this->MRImageProperty->SetLayerNumber(MRLayer);
  this->CTImageProperty = vtkSmartPointer<vtkImageProperty>::New();
  this->CTImageProperty->SetLayerNumber(CTLayer);
  this->MRVolumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
  this->ElectrodesProperty = vtkSmartPointer<vtkProperty>::New();

  // electrodes are rendered with a plain old actor
  this->ElectrodesActor = vtkSmartPointer<vtkActor>::New();
  vtkSmartPointer<vtkDataSetMapper> mapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  this->ElectrodesActor->SetMapper(mapper);
  this->ElectrodesActor->SetProperty(this->ElectrodesProperty);

  // volumes are cropped for efficient rendering
  this->BrainVolumeReslice = vtkSmartPointer<vtkImageReslice>::New();
  this->BrainVolumeReslice->SetInterpolationModeToCubic();
  this->HeadVolumeReslice = vtkSmartPointer<vtkImageReslice>::New();
  this->HeadVolumeReslice->SetInterpolationModeToCubic();

  // volume rendering is done with level-of-detail rendering
  vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> rayCastMapper;
  vtkSmartPointer<vtkVolumeTextureMapper3D> textureMapper;
  vtkSmartPointer<vtkPlane> clippingPlane =
    vtkSmartPointer<vtkPlane>::New();
  clippingPlane->SetNormal(1.0, 0.0, 0.0);
  clippingPlane->SetOrigin(0.0, 0.0, 0.0);
  int id;

  this->BrainVolume = vtkSmartPointer<vtkLODProp3D>::New();
  rayCastMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
  rayCastMapper->AutoAdjustSampleDistancesOff();
  rayCastMapper->LockSampleDistanceToInputSpacingOn();
  rayCastMapper->AddClippingPlane(clippingPlane);
  id = this->BrainVolume->AddLOD(rayCastMapper, this->MRVolumeProperty, 1.0);
  this->BrainVolume->DisableLOD(id);
  this->BrainVolumeLODIds.push_back(id);
  textureMapper = vtkSmartPointer<vtkVolumeTextureMapper3D>::New();
  textureMapper->AddClippingPlane(clippingPlane);
  id = this->BrainVolume->AddLOD(textureMapper, this->MRVolumeProperty, 0.1);
  this->BrainVolumeLODIds.push_back(id);

  this->HeadVolume = vtkSmartPointer<vtkLODProp3D>::New();
  rayCastMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
  rayCastMapper->AutoAdjustSampleDistancesOff();
  rayCastMapper->LockSampleDistanceToInputSpacingOn();
  id = this->HeadVolume->AddLOD(rayCastMapper, this->MRVolumeProperty, 1.0);
  this->BrainVolume->DisableLOD(id);
  this->HeadVolumeLODIds.push_back(id);
  this->HeadVolume->DisableLOD(id);
  textureMapper = vtkSmartPointer<vtkVolumeTextureMapper3D>::New();
  id = this->HeadVolume->AddLOD(textureMapper, this->MRVolumeProperty, 0.1);
  this->HeadVolumeLODIds.push_back(id);

  // generate an ortho slice for each orientation
  for (size_t i = 0; i < EpilepsyViewerDisplay::OrientationCount; ++i)
    {
    this->Slices.push_back(EpilepsyViewerDisplay::Slice(i));
    }

  // create the CT and MR layers
  this->AddLayer(this->CTImageProperty);
  this->AddLayer(this->MRImageProperty);

  // add everthing to the renderer
  this->MainRenderer->AddViewProp(this->BrainVolume);
  this->MainRenderer->AddViewProp(this->ElectrodesActor);
  size_t m = this->Slices.size();
  for (size_t j = 0; j < m; ++j)
    {
    this->MainRenderer->AddViewProp(this->Slices[j].Stack);
    }
}

//----------------------------------------------------------------------------
EpilepsyViewerDisplay::~EpilepsyViewerDisplay()
{
}

//----------------------------------------------------------------------------
void EpilepsyViewerDisplay::SetRenderWindow(vtkRenderWindow *renwin)
{
  this->RenderWindow = renwin;

  this->RenderWindow->AddRenderer(this->MainRenderer);

  // bind the OnRender() method to the window
  this->RenderWindow->AddObserver(
    vtkCommand::StartEvent, this, &EpilepsyViewerDisplay::OnRender);
}

//----------------------------------------------------------------------------
vtkRenderer *EpilepsyViewerDisplay::GetMainRenderer()
{
  return this->MainRenderer;
}

//----------------------------------------------------------------------------
vtkImageProperty *EpilepsyViewerDisplay::GetCTImageProperty()
{
  return this->CTImageProperty;
}

//----------------------------------------------------------------------------
vtkImageProperty *EpilepsyViewerDisplay::GetMRImageProperty()
{
  return this->MRImageProperty;
}

//----------------------------------------------------------------------------
vtkVolumeProperty *EpilepsyViewerDisplay::GetMRVolumeProperty()
{
  return this->MRVolumeProperty;
}

//----------------------------------------------------------------------------
vtkProperty *EpilepsyViewerDisplay::GetElectrodesProperty()
{
  return this->ElectrodesProperty;
}

//----------------------------------------------------------------------------
vtkImageSlice *EpilepsyViewerDisplay::GetImageSlice(int layer, int orientation)
{
  if (static_cast<size_t>(orientation) >= this->Slices.size())
    {
    return 0;
    }

  vtkImageStack *stack = this->Slices[orientation].Stack;
  vtkImageSliceCollection *images = stack->GetImages();

  vtkCollectionSimpleIterator pit;
  images->InitTraversal(pit);
  vtkImageSlice *image = 0;
  while ( (image = images->GetNextImage(pit)) != 0)
    {
    if (image->GetProperty()->GetLayerNumber() == layer)
      {
      break;
      }
    }

  return image;
}

//----------------------------------------------------------------------------
void EpilepsyViewerDisplay::SetData(EpilepsyViewerData *data)
{
  // the bounds of the brain are useful for many things
  vtkPolyData *mesh = data->GetMRBrainSurface();
  double bounds[6], center[4];
  mesh->GetBounds(bounds);
  center[0] = 0.5*(bounds[0] + bounds[1]);
  center[1] = 0.5*(bounds[2] + bounds[3]);
  center[2] = 0.5*(bounds[4] + bounds[5]);
  center[3] = 1.0;
  data->GetMRHeadMatrix()->MultiplyPoint(center, center);

  // the electrodes
  vtkDataSetMapper *mapper =
    vtkDataSetMapper::SafeDownCast(this->ElectrodesActor->GetMapper());
  mapper->SetInput(data->GetCTElectrodes());
  this->ElectrodesActor->SetUserMatrix(data->GetCTHeadMatrix());

  // the ortho planes
  vtkImageSlice *image = 0;
  int n = static_cast<int>(this->Slices.size());
  for (int i = 0; i < n; ++i)
    {
    double range[2];
    data->GetCTHeadAutoRange(range);
    image = this->GetImageSlice(CTLayer, i);
    image->GetMapper()->SetInput(data->GetCTHeadImage());
    image->GetProperty()->SetColorWindow(range[1] - range[0]);
    image->GetProperty()->SetColorLevel(0.5*(range[0] + range[1]));
    image->SetUserMatrix(data->GetCTHeadMatrix());

    data->GetMRHeadAutoRange(range);

    vtkSmartPointer<vtkTransform> resliceTransform =
      vtkSmartPointer<vtkTransform>::New();
    resliceTransform->PostMultiply();
    resliceTransform->Concatenate(data->GetMRHeadMatrix());
    resliceTransform->Inverse();
    resliceTransform->Concatenate(data->GetCTHeadMatrix());

    vtkSmartPointer<vtkImageReslice> mrReslice =
      vtkSmartPointer<vtkImageReslice>::New();
    mrReslice->SetInput(data->GetMRHeadImage());
    mrReslice->SetInformationInput(data->GetCTHeadImage());
    mrReslice->SetResliceTransform(resliceTransform);
    mrReslice->SetInterpolationModeToCubic();
    mrReslice->SetBackgroundLevel(range[0]);
    mrReslice->Update();

    image = this->GetImageSlice(MRLayer, i);
    image->GetMapper()->SetInput(data->GetCTHeadImage());
    image->GetMapper()->SetInput(mrReslice->GetOutput());
    image->GetProperty()->SetColorWindow(range[1] - range[0]);
    image->GetProperty()->SetColorLevel(0.5*(range[0] + range[1]));
    image->SetUserMatrix(data->GetCTHeadMatrix());
    }

  /*
  // crop the data before volume rendering
  vtkImageData *imageData = data->GetMRBrainImage();
  double origin[3], spacing[3];
  int extent[6];
  imageData->GetOrigin(origin);
  imageData->GetSpacing(spacing);
  imageData->GetWholeExtent(extent);
  for (int ii = 0; ii < 3; ++ii)
    {
    // add a 20 mm tolerance around the brain
    double minbound = bounds[2*ii] - 20.0;
    double maxbound = bounds[2*ii+1] + 20.0;
    int minext = static_cast<int>((minbound - origin[ii])/spacing[ii]);
    int maxext = static_cast<int>((maxbound - origin[ii])/spacing[ii]);
    extent[2*ii] = (minext < extent[2*ii] ? extent[2*ii] : minext);
    extent[2*ii+1] = (maxext > extent[2*ii+1] ? extent[2*ii+1] : maxext);
    }
  this->BrainVolumeReslice->SetInput(imageData);
  this->BrainVolumeReslice->SetOutputSpacing(spacing);
  this->BrainVolumeReslice->SetOutputOrigin(origin);
  this->BrainVolumeReslice->SetOutputExtent(extent);
  */

  // match the brain volume to the CT volume (temporary?)
  double range[2];
  data->GetMRHeadAutoRange(range);

  vtkSmartPointer<vtkTransform> resliceTransform =
    vtkSmartPointer<vtkTransform>::New();
  resliceTransform->PostMultiply();
  resliceTransform->Concatenate(data->GetMRHeadMatrix());
  resliceTransform->Inverse();
  resliceTransform->Concatenate(data->GetCTHeadMatrix());

  this->BrainVolumeReslice->SetInput(data->GetMRBrainImage());
  this->BrainVolumeReslice->SetInformationInput(data->GetCTHeadImage());
  this->BrainVolumeReslice->SetResliceTransform(resliceTransform);
  this->BrainVolumeReslice->SetInterpolationModeToCubic();
  this->BrainVolumeReslice->SetBackgroundLevel(range[0]);
  this->BrainVolumeReslice->Update();

  // the brain volume
  vtkLODProp3D *volume = this->BrainVolume;
  volume->SetUserMatrix(data->GetCTHeadMatrix());
  size_t m = this->HeadVolumeLODIds.size();
  for (size_t j = 0; j < m; ++j)
    {
    vtkAbstractVolumeMapper *vmapper;
    volume->GetLODMapper(this->BrainVolumeLODIds[j], &vmapper);
    if (vmapper)
      {
      vmapper->SetInputConnection(
        this->BrainVolumeReslice->GetOutputPort());
      }
    }

  // set up the volume property
  vtkSmartPointer<vtkColorTransferFunction> color =
    vtkSmartPointer<vtkColorTransferFunction>::New();
  vtkSmartPointer<vtkPiecewiseFunction> opacity =
    vtkSmartPointer<vtkPiecewiseFunction>::New();


  static double table[][5] = {
    { 0.00, 0.0, 0.0, 0.0, 0.0 },
    { 0.07, 0.4, 0.0, 0.1, 0.0 },
    { 0.48, 1.0, 0.7, 0.6, 0.2 },
    { 1.00, 1.0, 1.0, 0.9, 0.8 },
  };

  for (int i = 0; i < 4; ++i)
    {
    double x = table[i][0]*(range[1] - range[0]) + range[0];
    color->AddRGBPoint(x, table[i][1], table[i][2], table[i][3]);
    opacity->AddPoint(x, table[i][4]);
    }

  this->MRVolumeProperty->SetColor(color);
  this->MRVolumeProperty->SetScalarOpacity(opacity);
  this->MRVolumeProperty->SetInterpolationTypeToLinear();
  this->MRVolumeProperty->ShadeOff();

  // set the camera view using the brain surface
  this->MainRenderer->ResetCamera();
  vtkCamera *camera = this->MainRenderer->GetActiveCamera();
  camera->SetFocalPoint(center);
  center[0] -= 400.0;
  center[1] += 100.0;
  center[2] += 100.0;
  camera->SetPosition(center);
  camera->SetViewUp(0.0, 0.0, 1.0);
  camera->OrthogonalizeViewUp();
  camera->SetClippingRange(100, 700);
}

//----------------------------------------------------------------------------
void EpilepsyViewerDisplay::SetSliceOrientation(vtkMatrix4x4 *matrix)
{
  for (int i = 0; i < 3; ++i)
    {
    this->SetSliceOrientation(i, matrix);
    }
}

//----------------------------------------------------------------------------
void EpilepsyViewerDisplay::SetSliceOrientation(int i, vtkMatrix4x4 *matrix)
{
  // out-of-range check
  if (static_cast<size_t>(i) >= this->Slices.size())
    {
    return;
    }

  // get the plane to modify
  vtkPlane *plane = this->Slices[i].Plane;

  // get the orientation
  int o = (i < 3 ? i : 2);

  // homogeneous transformation of normal
  double normal[4] = { 0.0, 0.0, 0.0, 0.0 };
  double origin[3];
  plane->GetOrigin(origin);
  normal[i] = 1.0;
  normal[3] = -vtkMath::Dot(normal, origin);
  matrix->MultiplyPoint(normal, normal);
  plane->SetNormal(normal);
}


//----------------------------------------------------------------------------
void EpilepsyViewerDisplay::OnRender()
{
}

//----------------------------------------------------------------------------
void EpilepsyViewerDisplay::AddLayer(vtkImageProperty *imageProperty)
{
  vtkSmartPointer<vtkImageResliceMapper> mapper;
  vtkSmartPointer<vtkImageSlice> slice;

  size_t n = this->Slices.size();
  for (size_t i = 0; i < n; ++i)
    {
    slice = vtkSmartPointer<vtkImageSlice>::New();
    mapper = vtkSmartPointer<vtkImageResliceMapper>::New();
    mapper->SetSlicePlane(this->Slices[i].Plane);
    slice->SetMapper(mapper);
    slice->SetProperty(imageProperty);
    this->Slices[i].Stack->AddImage(slice);
    }
}
