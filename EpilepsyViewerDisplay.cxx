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
#include <vtkPlaneCollection.h>
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
  this->ElectrodesProperty->SetColor(0.4, 0.7, 1.0);

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
  //this->BrainVolume->DisableLOD(id);
  this->BrainVolumeLODIds.push_back(id);

  this->HeadVolume = vtkSmartPointer<vtkLODProp3D>::New();
  rayCastMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
  rayCastMapper->AutoAdjustSampleDistancesOff();
  rayCastMapper->LockSampleDistanceToInputSpacingOn();
  id = this->HeadVolume->AddLOD(rayCastMapper, this->MRVolumeProperty, 1.0);
  this->HeadVolume->DisableLOD(id);
  this->HeadVolumeLODIds.push_back(id);
  textureMapper = vtkSmartPointer<vtkVolumeTextureMapper3D>::New();
  id = this->HeadVolume->AddLOD(textureMapper, this->MRVolumeProperty, 0.1);
  //this->HeadVolume->DisableLOD(id);
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
  // get the center of the brain mesh, us it as the focus
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

  // transformation from CT to MR
  vtkSmartPointer<vtkTransform> transformCTtoMR =
    vtkSmartPointer<vtkTransform>::New();
  transformCTtoMR->PostMultiply();
  transformCTtoMR->Concatenate(data->GetMRHeadMatrix());
  transformCTtoMR->Inverse();
  transformCTtoMR->Concatenate(data->GetCTHeadMatrix());

  // get CT bounds in MR coords, be careful about half-voxel border
  double spacing[3], origin[3];
  double ctbounds[6];
  data->GetCTHeadImage()->GetBounds(ctbounds);
  data->GetCTHeadImage()->GetSpacing(spacing);
  for (int i = 0; i < 3; ++i)
    {
    double b = 0.5*fabs(spacing[i]);
    ctbounds[2*i] -= b;
    ctbounds[2*i + 1] += b;
    }
  EpilepsyViewerDisplay::TransformBounds(
    transformCTtoMR->GetMatrix(), ctbounds, bounds);

  // pad the MR to avoid showing CT past its bounds
  // (eventually ImageResliceMapper should have a setting that
  // will make this unnecessary).
  int extent[6];
  data->GetMRHeadImage()->GetSpacing(spacing);
  data->GetMRHeadImage()->GetOrigin(origin);
  data->GetMRHeadImage()->GetWholeExtent(extent);

  for (int i = 0; i < 3; ++i)
    {
    if (spacing[i] < 0)
      {
      origin[i] += spacing[i]*(extent[2*i+1] - extent[2*i]);
      spacing[i] = -spacing[i];
      }

    double b = 0.5*spacing[i];
    bounds[2*i] += b;
    bounds[2*i + 1] -= b;

    origin[i] += vtkMath::Floor((bounds[2*i]-origin[i])/spacing[i])*spacing[i];
    extent[2*i] = 0;
    extent[2*i + 1] = vtkMath::Ceil((bounds[2*i + 1] - origin[i])/spacing[i]);
    }

/*
  vtkSmartPointer<vtkImageReslice> padFilter =
    vtkSmartPointer<vtkImageReslice>::New();
  padFilter->SetInput(data->GetMRHeadImage());
  padFilter->SetBackgroundLevel(range[0]);
  padFilter->SetOutputSpacing(spacing);
  padFilter->SetOutputOrigin(origin);
  padFilter->SetOutputExtent(extent);
  padFilter->Update();
*/
/*
  // also pad the volume, for the sake of its clipping cube
  this->BrainVolumeReslice->SetInput(data->GetMRBrainImage());
  this->BrainVolumeReslice->SetBackgroundLevel(range[0]);
  this->BrainVolumeReslice->SetOutputSpacing(spacing);
  this->BrainVolumeReslice->SetOutputOrigin(origin);
  this->BrainVolumeReslice->SetOutputExtent(extent);
  this->BrainVolumeReslice->Update();
*/
  // the ortho planes
  vtkImageSlice *image = 0;
  int n = static_cast<int>(this->Slices.size());
  for (int i = 0; i < n; ++i)
    {
    // set the plane positions
    this->Slices[i].Plane->SetOrigin(center);

    // set up the planes for the CT
    double range[2];
    data->GetCTHeadAutoRange(range);
    image = this->GetImageSlice(CTLayer, i);
    image->GetMapper()->SetInput(data->GetCTHeadImage());
    image->GetMapper()->BorderOn();
    image->GetProperty()->SetColorWindow(range[1] - range[0]);
    image->GetProperty()->SetColorLevel(0.5*(range[0] + range[1]));
    image->SetUserMatrix(data->GetCTHeadMatrix());

    // use CT to generate clipping planes
    vtkSmartPointer<vtkPlaneCollection> clippingPlanes =
      vtkSmartPointer<vtkPlaneCollection>::New();
    EpilepsyViewerDisplay::GenerateClippingPlanes(
      clippingPlanes, data->GetCTHeadMatrix(), ctbounds);

    // get auto-range for MR
    data->GetMRHeadAutoRange(range);

    image = this->GetImageSlice(MRLayer, i);
    image->GetMapper()->SetClippingPlanes(clippingPlanes);
    //image->GetMapper()->SetInput(padFilter->GetOutput());
    image->GetMapper()->SetInput(data->GetMRHeadImage());
    image->GetMapper()->BorderOn();
    image->GetProperty()->SetColorWindow(range[1] - range[0]);
    image->GetProperty()->SetColorLevel(0.5*(range[0] + range[1]));
    image->SetUserMatrix(data->GetMRHeadMatrix());
    }

  // the brain volume
  vtkLODProp3D *volume = this->BrainVolume;
  volume->SetUserMatrix(data->GetMRHeadMatrix());

  // use CT to generate clipping planes
  vtkSmartPointer<vtkPlaneCollection> clippingPlanes =
    vtkSmartPointer<vtkPlaneCollection>::New();
  EpilepsyViewerDisplay::GenerateClippingPlanes(
    clippingPlanes, data->GetCTHeadMatrix(), ctbounds);
  clippingPlanes->GetItem(5)->SetOrigin(center[0], center[1], center[2] + 30);

  size_t m = this->HeadVolumeLODIds.size();
  for (size_t j = 0; j < m; ++j)
    {
    vtkAbstractVolumeMapper *vmapper;
    volume->GetLODMapper(this->BrainVolumeLODIds[j], &vmapper);
    if (vmapper)
      {
      vmapper->SetClippingPlanes(clippingPlanes);
      //vmapper->SetInputConnection(
      //  this->BrainVolumeReslice->GetOutputPort());
      vmapper->SetInput(data->GetMRBrainImage());
      }
    }

  // set up the volume property
  vtkSmartPointer<vtkColorTransferFunction> color =
    vtkSmartPointer<vtkColorTransferFunction>::New();
  vtkSmartPointer<vtkPiecewiseFunction> opacity =
    vtkSmartPointer<vtkPiecewiseFunction>::New();

  double range[2];
  data->GetMRHeadAutoRange(range);

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
  center[0] -= 600.0;
  center[1] += 150.0;
  center[2] += 150.0;
  camera->SetPosition(center);
  camera->SetViewUp(0.0, 0.0, 1.0);
  camera->OrthogonalizeViewUp();
  camera->SetClippingRange(200, 1000);
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

//----------------------------------------------------------------------------
void EpilepsyViewerDisplay::GenerateClippingPlanes(
  vtkPlaneCollection *planes, vtkMatrix4x4 *matrix, const double bounds[6])
{
  double *pointMatrix = *matrix->Element;
  double normalMatrix[16];
  vtkMatrix4x4::Invert(pointMatrix, normalMatrix);
  vtkMatrix4x4::Transpose(normalMatrix, normalMatrix);

  double center[3];
  center[0] = 0.5*(bounds[0] + bounds[1]);
  center[1] = 0.5*(bounds[2] + bounds[3]);
  center[2] = 0.5*(bounds[4] + bounds[5]);

  planes->RemoveAllItems();

  for (int i = 0; i < 6; ++i)
    {
    int j = (i >> 1);

    double origin[4], normal[4];
    origin[0] = center[0];
    origin[1] = center[1];
    origin[2] = center[2];
    origin[3] = 1.0;

    origin[j] = bounds[i];

    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    normal[3] = 0.0;

    normal[j] = 1 - ((i & 1) << 1);

    normal[3] = -normal[j]*origin[j];

    vtkMatrix4x4::MultiplyPoint(pointMatrix, origin, origin);
    origin[0] /= origin[3];
    origin[1] /= origin[3];
    origin[2] /= origin[3];

    vtkMatrix4x4::MultiplyPoint(normalMatrix, normal, normal);
    double norm =
      sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    normal[0] /= norm;
    normal[1] /= norm;
    normal[2] /= norm;

    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(origin);
    plane->SetNormal(normal);

    planes->AddItem(plane);
    }
}

//----------------------------------------------------------------------------
void EpilepsyViewerDisplay::TransformBounds(
  vtkMatrix4x4 *matrix, const double inbounds[6], double outbounds[6])
{
  double xmin = VTK_DOUBLE_MAX;
  double xmax = VTK_DOUBLE_MIN;
  double ymin = VTK_DOUBLE_MAX;
  double ymax = VTK_DOUBLE_MIN;
  double zmin = VTK_DOUBLE_MAX;
  double zmax = VTK_DOUBLE_MIN;

  for (int n = 0; n < 8; ++n)
    {
    double corner[4];
    corner[0] = inbounds[(n & 1)];
    corner[1] = inbounds[2 + ((n & 2) >> 1)];
    corner[2] = inbounds[4 + ((n & 4) >> 2)];
    corner[3] = 1.0;

    matrix->MultiplyPoint(corner, corner);

    double x = corner[0] / corner[3];
    xmin = (x > xmin ? xmin : x);
    xmax = (x < xmax ? xmax : x);
    double y = corner[1] / corner[3];
    ymin = (y > ymin ? ymin : y);
    ymax = (y < ymax ? ymax : y);
    double z = corner[2] / corner[3];
    zmin = (z > zmin ? zmin : z);
    zmax = (z < zmax ? zmax : z);
    }

  outbounds[0] = xmin;
  outbounds[1] = xmax;
  outbounds[2] = ymin;
  outbounds[3] = ymax;
  outbounds[4] = zmin;
  outbounds[5] = zmax;
}
