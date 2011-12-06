/*=========================================================================

Program:   EpilepsyViewerData
Module:    EpilepsyViewerData.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

#include "EpilepsyViewerData.h"

// from AIRS
#include <vtkImageMRIBrainExtractor.h>
#include <vtkImageRegistration.h>
#include <vtkImageFastBlur.h>

// from VTK
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkImageReslice.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMarchingCubes.h>
#include <vtkROIStencilSource.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkStripper.h>
#include <vtkTimerLog.h>
#include <vtkImageSincInterpolator.h>
#include <vtkMath.h>

#include <vtkMINCImageReader.h>
#include <vtkDICOMImageReader.h>

#include <vtksys/SystemTools.hxx>
#include <vtksys/Glob.hxx>

// from the standard library
#include <stddef.h>
#include <stdarg.h>

//----------------------------------------------------------------------------
EpilepsyViewerData::EpilepsyViewerData()
{
  this->CTHeadImage = vtkSmartPointer<vtkImageData>::New();
  this->MRHeadImage = vtkSmartPointer<vtkImageData>::New();
  this->MRBrainImage = vtkSmartPointer<vtkImageData>::New();

  this->CTHeadMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  this->MRHeadMatrix = vtkSmartPointer<vtkMatrix4x4>::New();

  this->MRBrainSurface = vtkSmartPointer<vtkPolyData>::New();
  this->CTElectrodes = vtkSmartPointer<vtkPolyData>::New();

  this->Timer = vtkSmartPointer<vtkTimerLog>::New();
}

//----------------------------------------------------------------------------
EpilepsyViewerData::~EpilepsyViewerData()
{
}

//----------------------------------------------------------------------------
void EpilepsyViewerData::ReportError(const char *format, ...)
{
  char errmessage[1024];
  va_list ap;
  va_start(ap, format);
  int nc = vsnprintf(errmessage, 1024, format, ap);
  va_end(ap);

  if (nc <= 0)
    {
    cerr << "[error message format error]\n";
    }
  else
    {
    cerr << errmessage << "\n";
    }
}

//----------------------------------------------------------------------------
vtkImageData *EpilepsyViewerData::GetCTHeadImage()
{
  return this->CTHeadImage;
}

//----------------------------------------------------------------------------
vtkMatrix4x4 *EpilepsyViewerData::GetCTHeadMatrix()
{
  return this->CTHeadMatrix;
}

//----------------------------------------------------------------------------
vtkImageData *EpilepsyViewerData::GetMRHeadImage()
{
  return this->MRHeadImage;
}

//----------------------------------------------------------------------------
vtkImageData *EpilepsyViewerData::GetMRBrainImage()
{
  return this->MRBrainImage;
}

//----------------------------------------------------------------------------
vtkMatrix4x4 *EpilepsyViewerData::GetMRHeadMatrix()
{
  return this->MRHeadMatrix;
}

//----------------------------------------------------------------------------
vtkPolyData *EpilepsyViewerData::GetCTElectrodes()
{
  return this->CTElectrodes;
}

//----------------------------------------------------------------------------
vtkPolyData *EpilepsyViewerData::GetMRBrainSurface()
{
  return this->MRBrainSurface;
}

//----------------------------------------------------------------------------
bool EpilepsyViewerData::LoadFromDataDirectory(const char *dirname)
{
  if (!vtksys::SystemTools::FileExists(dirname))
    {
    this->ReportError("Directory does not exist: %s", dirname);
    return false;
    }

  if (!vtksys::SystemTools::FileIsDirectory(dirname))
    {
    this->ReportError("Not a directory: %s", dirname);
    return false;
    }

  std::vector<std::string> components;
  vtksys::SystemTools::SplitPath(dirname, components);

  components.push_back("*CT.mnc");
  std::string mincCT = vtksys::SystemTools::JoinPath(components);
  components.pop_back();

  components.push_back("*MR.mnc");
  std::string mincMR = vtksys::SystemTools::JoinPath(components);
  components.pop_back();

  vtksys::Glob glob;
  if (glob.FindFiles(mincCT))
    {
    std::vector<std::string> &files = glob.GetFiles();
    if (files.size() == 1)
      {
      mincCT = files.front();
      }
    }

  if (glob.FindFiles(mincMR))
    {
    std::vector<std::string> &files = glob.GetFiles();
    if (files.size() == 1)
      {
      mincMR = glob.GetFiles().front();
      }
    }

  if (vtksys::SystemTools::FileExists(mincCT.c_str(), true))
    {
    bool success = this->ReadMINCImage(
      this->CTHeadImage, this->CTHeadMatrix,
      this->CTHeadFullRange, this->CTHeadAutoRange,
      mincCT.c_str());

    if (!success)
      {
      this->ReportError("Bad file: %s", mincCT.c_str());
      return false;
      }
    }

  if (vtksys::SystemTools::FileExists(mincMR.c_str(), true))
    {
    bool success = this->ReadMINCImage(
      this->MRHeadImage, this->MRHeadMatrix,
      this->MRHeadFullRange, this->MRHeadAutoRange,
      mincMR.c_str());

    if (!success)
      {
      this->ReportError("Bad file: %s", mincMR.c_str());
      return false;
      }
    }

  return true;
}

//----------------------------------------------------------------------------
bool EpilepsyViewerData::RegisterMRHeadToCTHead()
{
  // -------------------------------------------------------
  // parameters for registration

  int interpolatorType = vtkImageRegistration::Linear;
  double transformTolerance = 0.1; // tolerance on transformation result
  int numberOfBins = 64; // for Mattes' mutual information
  double initialBlurFactor = 4.0;

  // -------------------------------------------------------
  // the data to be registered

  vtkImageData *targetImage = this->CTHeadImage;
  vtkImageData *sourceImage = this->MRHeadImage;
  vtkMatrix4x4 *targetMatrix = this->CTHeadMatrix;
  vtkMatrix4x4 *sourceMatrix = this->MRHeadMatrix;

  // -------------------------------------------------------
  // prepare for registration

  // get information about the images
  double targetSpacing[3], sourceSpacing[3];
  targetImage->GetSpacing(targetSpacing);
  sourceImage->GetSpacing(sourceSpacing);

  for (int jj = 0; jj < 3; jj++)
    {
    targetSpacing[jj] = fabs(targetSpacing[jj]);
    sourceSpacing[jj] = fabs(sourceSpacing[jj]);
    }

  double minSpacing = sourceSpacing[0];
  if (minSpacing > sourceSpacing[1])
    {
    minSpacing = sourceSpacing[1];
    }
  if (minSpacing > sourceSpacing[2])
    {
    minSpacing = sourceSpacing[2];
    }

  // blur source image with Hamming-windowed sinc
  vtkSmartPointer<vtkImageSincInterpolator> sourceBlurKernel =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  sourceBlurKernel->SetWindowFunctionToHamming();

  // reduce the source resolution
  vtkSmartPointer<vtkImageFastBlur> sourceBlur =
    vtkSmartPointer<vtkImageFastBlur>::New();
  sourceBlur->SetInput(sourceImage);
  sourceBlur->SetResizeMethodToOutputSpacing();
  sourceBlur->SetInterpolator(sourceBlurKernel);
  sourceBlur->InterpolateOn();

  // blur target with Hamming-windowed sinc
  vtkSmartPointer<vtkImageSincInterpolator> targetBlurKernel =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  targetBlurKernel->SetWindowFunctionToHamming();

  // keep target at full resolution
  vtkSmartPointer<vtkImageFastBlur> targetBlur =
    vtkSmartPointer<vtkImageFastBlur>::New();
  targetBlur->SetInput(targetImage);
  targetBlur->SetResizeMethodToOutputSpacing();
  targetBlur->SetInterpolator(targetBlurKernel);
  targetBlur->InterpolateOn();

  // get the initial transformation
  vtkSmartPointer<vtkMatrix4x4> matrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->DeepCopy(targetMatrix);
  matrix->Invert();
  vtkMatrix4x4::Multiply4x4(matrix, sourceMatrix, matrix);

  // set up the registration
  vtkSmartPointer<vtkImageRegistration> registration =
    vtkSmartPointer<vtkImageRegistration>::New();
  registration->SetTargetImageInputConnection(targetBlur->GetOutputPort());
  registration->SetSourceImageInputConnection(sourceBlur->GetOutputPort());
  registration->SetInitializerTypeToCentered();
  registration->SetTransformTypeToRigid();
  registration->SetMetricTypeToNormalizedMutualInformation();
  registration->SetInterpolatorType(interpolatorType);
  registration->SetJointHistogramSize(numberOfBins,numberOfBins);
  registration->SetMetricTolerance(1e-4);
  registration->SetTransformTolerance(transformTolerance);
  registration->SetMaximumNumberOfIterations(500);

  // -------------------------------------------------------
  // start the timer
  double startTime = this->Timer->GetUniversalTime();
  double lastTime = startTime;

  // -------------------------------------------------------
  // do the registration

  // the registration starts at low-resolution
  double blurFactor = initialBlurFactor;
  // two stages for each resolution:
  // first without interpolation, and then with interpolation
  int stage = 0;
  // will be set to "true" when registration is initialized
  bool initialized = false;

  for (;;)
    {
    if (stage == 0)
      {
      registration->SetInterpolatorTypeToNearest();
      registration->SetTransformTolerance(minSpacing*blurFactor);
      }
    else
      {
      registration->SetInterpolatorType(interpolatorType);
      registration->SetTransformTolerance(transformTolerance*blurFactor);
      }
    if (blurFactor < 1.1)
      {
      // full resolution: no blurring or resampling
      sourceBlur->InterpolateOff();
      sourceBlur->SetOutputSpacing(sourceSpacing);
      sourceBlur->Update();

      targetBlur->InterpolateOff();
      targetBlur->SetOutputSpacing(targetSpacing);
      targetBlur->Update();
      }
    else
      {
      // reduced resolution: set the blurring
      double spacing[3];
      for (int j = 0; j < 3; j++)
        {
        spacing[j] = blurFactor*minSpacing;
        if (spacing[j] < sourceSpacing[j])
          {
          spacing[j] = sourceSpacing[j];
          }
        }

      sourceBlurKernel->SetBlurFactors(
        spacing[0]/sourceSpacing[0],
        spacing[1]/sourceSpacing[1],
        spacing[2]/sourceSpacing[2]);

      sourceBlur->SetOutputSpacing(spacing);
      sourceBlur->Update();

      targetBlurKernel->SetBlurFactors(
        blurFactor*minSpacing/targetSpacing[0],
        blurFactor*minSpacing/targetSpacing[1],
        blurFactor*minSpacing/targetSpacing[2]);

      targetBlur->Update();

      if (stage == 0)
        {
        double nTime = this->Timer->GetUniversalTime();
        cout << "blur " << blurFactor << " blurring took "
             << (nTime - lastTime) << " seconds" << endl;
        lastTime = nTime;
        }
      }

    if (initialized)
      {
      // re-initialize with the matrix from the previous step
      registration->SetInitializerTypeToNone();
      matrix->DeepCopy(registration->GetTransform()->GetMatrix());
      }

    registration->Initialize(matrix);

    initialized = true;

    while (registration->Iterate())
      {
      // will iterate until convergence or failure
      vtkMatrix4x4::Multiply4x4(
        targetMatrix,registration->GetTransform()->GetMatrix(),sourceMatrix);
      sourceMatrix->Modified();
      }

    double newTime = this->Timer->GetUniversalTime();
    cout << "blur " << blurFactor << " stage " << stage << " took "
         << (newTime - lastTime) << "s and "
         << registration->GetNumberOfEvaluations() << " evaluations" << endl;
    lastTime = newTime;

    // check to see how much the registration changed
    matrix->Invert();
    vtkMatrix4x4::Multiply4x4(
      matrix, registration->GetTransform()->GetMatrix(), matrix);
    double dt = sqrt(matrix->Element[0][3]*matrix->Element[0][3] +
                     matrix->Element[1][3]*matrix->Element[1][3] +
                     matrix->Element[2][3]*matrix->Element[2][3]);
    double dr = (matrix->Element[0][0] +
                 matrix->Element[1][1] +
                 matrix->Element[2][2])/3;
    dr = ((dr > 0) ? dr : 0);
    dr = ((dr < 1) ? dr : 1);
    dr = acos(sqrt(dr))*180/vtkMath::DoublePi();
    cout << "blur " << blurFactor << " stage " << stage
         << " changed result by " << dt << " mm " << dr << " degrees" << endl;

    // prepare for next iteration
    if (stage == 1)
      {
      blurFactor /= 2.0;
      if (blurFactor < 0.9)
        {
        break;
        }
      }
    stage = (stage + 1) % 2;
    }

  cout << "registration took " << (lastTime - startTime) << "s" << endl;

  return true;
}

//----------------------------------------------------------------------------
bool EpilepsyViewerData::ExtractMRBrain()
{
  // start the timer
  double startTime = this->Timer->GetUniversalTime();

  vtkSmartPointer<vtkImageMRIBrainExtractor> extractor =
    vtkSmartPointer<vtkImageMRIBrainExtractor>::New();

  double bt = 0.70;

  extractor->SetInput(this->MRHeadImage);
  extractor->SetRMin(8.0);
  extractor->SetRMax(10.0);
  extractor->SetD1(7.0);
  extractor->SetD2(3.0);
  extractor->SetBT(bt);
  extractor->Update();
  vtkImageData *image = extractor->GetOutput();
  vtkPolyData *mesh = extractor->GetBrainMesh();

  // get the data
  this->MRBrainImage->CopyStructure(image);
  this->MRBrainImage->GetPointData()->PassData(image->GetPointData());

  // get the mesh
  this->MRBrainSurface->SetPoints(mesh->GetPoints());
  this->MRBrainSurface->SetPolys(mesh->GetPolys());

  double endTime = this->Timer->GetUniversalTime();
  cerr << "brain extraction took " << (endTime - startTime)
       << " seconds" << endl;

  return true;
}

//----------------------------------------------------------------------------
bool EpilepsyViewerData::ExtractCTElectrodes()
{
  // start the timer
  double startTime = this->Timer->GetUniversalTime();

  // the image to extract the electrodes from
  vtkImageData *image = this->CTHeadImage;

  // get the range of the data, the electrodes will be above
  // the autoRange because they make up a very small percentage
  // of the voxels in the image
  vtkSmartPointer<vtkImageHistogramStatistics> rangeFinder =
    vtkSmartPointer<vtkImageHistogramStatistics>::New();

  rangeFinder->SetInput(this->CTHeadImage);
  rangeFinder->Update();
  double autoRange[2];
  rangeFinder->GetAutoRange(autoRange);
  double skullValue = autoRange[1];
  double maxValue = rangeFinder->GetMaximum();

  // set the isovalue to be higher the skull, lower than electrodes
  double value = (0.5*skullValue + 0.5*maxValue);

  // get the image spacing, convert to isotropic
  double spacing[3];
  double origin[3];
  int extent[6];
  image->GetSpacing(spacing);
  image->GetOrigin(origin);
  image->GetWholeExtent(extent);
  double zsize = (extent[5] - extent[4])*spacing[2];
  double tol = zsize*1e-3;
  spacing[2] = spacing[0];
  extent[5] = static_cast<int>(zsize/spacing[2] + tol);

  // find the center of the brain mesh
  double bounds[6];
  double center[3];
  this->MRBrainSurface->GetBounds(bounds);
  center[0] = 0.5*(bounds[0] + bounds[1]);
  center[1] = 0.5*(bounds[2] + bounds[3]);
  center[2] = 0.5*(bounds[4] + bounds[5]);

  // transform to put the brain mesh into the CT data space
  vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();
  vtkSmartPointer<vtkMatrix4x4> matrix =
    vtkSmartPointer<vtkMatrix4x4>::New();

  vtkMatrix4x4::Invert(this->CTHeadMatrix, matrix);
  matrix->Multiply4x4(this->MRHeadMatrix, matrix, matrix);
  transform->Concatenate(matrix);

  // expand by a 5% tolerance to get all electrodes
  double f = 1.05;
  transform->TransformPoint(center, center);
  transform->PostMultiply();
  transform->Translate(-center[0], -center[1], -center[2]);
  transform->Scale(f, f, f);
  transform->Translate(center[0], center[1], center[2]);

  // apply the transform to the mesh
  vtkSmartPointer<vtkTransformPolyDataFilter> tfilter =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();

  tfilter->SetInput(this->MRBrainSurface);
  tfilter->SetTransform(transform);

  // convert the brain mesh into a stencil
  vtkSmartPointer<vtkPolyDataToImageStencil> makeStencil =
    vtkSmartPointer<vtkPolyDataToImageStencil>::New();

  makeStencil->SetInputConnection(tfilter->GetOutputPort());
  makeStencil->SetOutputSpacing(spacing);
  makeStencil->SetOutputOrigin(origin);
  makeStencil->SetOutputWholeExtent(extent);
  makeStencil->Update();

  // apply the stencil while resampling the image
  vtkSmartPointer<vtkImageSincInterpolator> resampleKernel =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  vtkSmartPointer<vtkImageReslice> stencil =
    vtkSmartPointer<vtkImageReslice>::New();

  // use high-quality filter when creating isotropic spacing
  resampleKernel->SetWindowFunctionToBlackman();

  stencil->SetInput(image);
  stencil->SetInterpolator(resampleKernel);
  stencil->SetStencil(makeStencil->GetOutput());
  stencil->SetOutputSpacing(spacing);
  stencil->SetOutputOrigin(origin);
  stencil->SetOutputExtent(extent);

  // run marching cubes
  vtkSmartPointer<vtkMarchingCubes> cubes =
    vtkSmartPointer<vtkMarchingCubes>::New();

  cubes->SetInputConnection(stencil->GetOutputPort());
  cubes->SetValue(0, value);
  cubes->ComputeNormalsOn();

  // render faster by stripping
  vtkSmartPointer<vtkStripper> stripper =
    vtkSmartPointer<vtkStripper>::New();

  stripper->SetInputConnection(cubes->GetOutputPort());

  stripper->Update();
  vtkPolyData *mesh = stripper->GetOutput();

  // get the mesh, ignore the scalars
  this->CTElectrodes->SetPoints(mesh->GetPoints());
  this->CTElectrodes->SetPolys(mesh->GetPolys());
  this->CTElectrodes->SetStrips(mesh->GetStrips());
  this->CTElectrodes->GetPointData()->SetNormals(
    mesh->GetPointData()->GetNormals());

  double endTime = this->Timer->GetUniversalTime();
  cerr << "electrode extraction took " << (endTime - startTime)
       << " seconds" << endl;

  return true;
}

//----------------------------------------------------------------------------
bool EpilepsyViewerData::ReadDICOMSeries(
  vtkImageData *data, vtkMatrix4x4 *matrix,
  double fullRange[2], double autoRange[2],
  const char *directoryName)
{
  // read the image
  vtkSmartPointer<vtkDICOMImageReader> reader =
    vtkSmartPointer<vtkDICOMImageReader>::New();

  reader->SetDirectoryName(directoryName);
  reader->Update();
  vtkImageData *image = reader->GetOutput();

  // compute the new image info for DICOM-style ordering
  double spacing[3];
  double origin[3];
  int extent[6];
  double bounds[6];
  image->GetSpacing(spacing);
  image->GetOrigin(origin);
  image->GetWholeExtent(extent);
  for (int i = 0; i < 3; ++i)
    {
    double b1 = extent[2*i]*spacing[i] + origin[i];
    double b2 = extent[2*i+1]*spacing[i] + origin[i];
    b1 = (i == 0 ? b1 : -b1); // flip if not X
    b2 = (i == 0 ? b2 : -b2); // flip if not X
    bounds[2*i] = (b1 < b2 ? b1 : b2);
    bounds[2*i+1] = (b1 < b2 ? b2 : b1);
    spacing[i] = fabs(spacing[i]);
    origin[i] = bounds[2*i];
    // reduce bounds by 2% in X and Y for use in cylinder generation
    double bl = (i == 2 ? 0.0 : 0.01*(bounds[2*i+1] - bounds[2*i]));
    bounds[2*i] += bl;
    bounds[2*i+1] -= bl;
    }

  // extract just the reconstructed portion of CT image
  vtkSmartPointer<vtkROIStencilSource> cylinder =
    vtkSmartPointer<vtkROIStencilSource>::New();

  cylinder->SetShapeToCylinderZ();
  cylinder->SetOutputSpacing(spacing);
  cylinder->SetOutputOrigin(origin);
  cylinder->SetOutputWholeExtent(extent);
  cylinder->SetBounds(bounds);

  // get the range within the cylinder
  vtkSmartPointer<vtkImageHistogramStatistics> rangeFinder =
    vtkSmartPointer<vtkImageHistogramStatistics>::New();

  rangeFinder->SetInputConnection(reader->GetOutputPort());
  rangeFinder->SetStencil(cylinder->GetOutput());
  rangeFinder->Update();

  // if the data range in the cylinder is significantly smaller than
  // the data range of entire image, then image has a circular ROI
  image->GetScalarRange(fullRange);
  rangeFinder->GetAutoRange(autoRange);
  double newRange[2];
  newRange[0] = rangeFinder->GetMinimum();
  newRange[1] = rangeFinder->GetMaximum();
  bool hasCircularROI =
    ((fullRange[1] - fullRange[0] + 1)/(newRange[1] - newRange[0] + 1) > 1.1);

  // the reader flips the image and reverses the ordering, so undo these
  vtkSmartPointer<vtkImageReslice> flip =
    vtkSmartPointer<vtkImageReslice>::New();

  if (hasCircularROI)
    {
    fullRange[0] = newRange[0];
    fullRange[1] = newRange[1];
    flip->SetStencil(cylinder->GetOutput());
    std::vector<std::string> components;
    vtksys::SystemTools::SplitPath(directoryName, components);
    cout << "extracting circular ROI from " << components.back() << endl;
    }
  flip->SetInputConnection(reader->GetOutputPort());
  flip->SetResliceAxesDirectionCosines(1,0,0, 0,-1,0, 0,0,-1);
  flip->SetOutputOrigin(origin);
  flip->SetOutputExtent(extent);
  flip->Update();

  image = flip->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());
  data->SetOrigin(0,0,0);

  // generate the matrix
  float *position = reader->GetImagePositionPatient();
  float *orientation = reader->GetImageOrientationPatient();
  float *xdir = &orientation[0];
  float *ydir = &orientation[3];
  float zdir[3];
  vtkMath::Cross(xdir, ydir, zdir);

  for (int i = 0; i < 3; i++)
    {
    matrix->Element[i][0] = xdir[i];
    matrix->Element[i][1] = ydir[i];
    matrix->Element[i][2] = zdir[i];
    matrix->Element[i][3] = position[i];
    }
  matrix->Element[3][0] = 0;
  matrix->Element[3][1] = 0;
  matrix->Element[3][2] = 0;
  matrix->Element[3][3] = 1;
  matrix->Modified();

  return true;
}

//----------------------------------------------------------------------------
bool EpilepsyViewerData::ReadMINCImage(
  vtkImageData *data, vtkMatrix4x4 *matrix,
  double fullRange[2], double autoRange[2], const char *fileName)
{
  // read the image
  vtkSmartPointer<vtkMINCImageReader> reader =
    vtkSmartPointer<vtkMINCImageReader>::New();

  reader->SetFileName(fileName);
  reader->Update();
  vtkImageData *image = reader->GetOutput();

  // compute the new image info for DICOM-style ordering
  double spacing[3];
  double origin[3];
  int extent[6];
  double bounds[6];
  image->GetSpacing(spacing);
  image->GetOrigin(origin);
  image->GetWholeExtent(extent);
  for (int i = 0; i < 3; ++i)
    {
    double b1 = extent[2*i]*spacing[i] + origin[i];
    double b2 = extent[2*i+1]*spacing[i] + origin[i];
    b1 = (i == 2 ? b1 : -b1); // flip if not Z
    b2 = (i == 2 ? b2 : -b2); // flip if not Z
    bounds[2*i] = (b1 < b2 ? b1 : b2);
    bounds[2*i+1] = (b1 < b2 ? b2 : b1);
    spacing[i] = fabs(spacing[i]);
    origin[i] = bounds[2*i];
    // reduce bounds by 2% in X and Y for use in cylinder generation
    double bl = (i == 2 ? 0.0 : 0.01*(bounds[2*i+1] - bounds[2*i]));
    bounds[2*i] += bl;
    bounds[2*i+1] -= bl;
    }

  // extract just the reconstructed portion of CT image
  vtkSmartPointer<vtkROIStencilSource> cylinder =
    vtkSmartPointer<vtkROIStencilSource>::New();

  cylinder->SetShapeToCylinderZ();
  cylinder->SetOutputSpacing(spacing);
  cylinder->SetOutputOrigin(origin);
  cylinder->SetOutputWholeExtent(extent);
  cylinder->SetBounds(bounds);

  // get the range within the cylinder
  vtkSmartPointer<vtkImageHistogramStatistics> rangeFinder =
    vtkSmartPointer<vtkImageHistogramStatistics>::New();

  rangeFinder->SetInputConnection(reader->GetOutputPort());
  rangeFinder->SetStencil(cylinder->GetOutput());
  rangeFinder->Update();

  // if the data range in the cylinder is significantly smaller than
  // the data range of entire image, then image has a circular ROI
  image->GetScalarRange(fullRange);
  rangeFinder->GetAutoRange(autoRange);
  double newRange[2];
  newRange[0] = rangeFinder->GetMinimum();
  newRange[1] = rangeFinder->GetMaximum();
  bool hasCircularROI =
    ((fullRange[1] - fullRange[0] + 1)/(newRange[1] - newRange[0] + 1) > 1.1);

  // flip the image rows into a DICOM-style ordering
  vtkSmartPointer<vtkImageReslice> flip =
    vtkSmartPointer<vtkImageReslice>::New();

  if (hasCircularROI)
    {
    fullRange[0] = newRange[0];
    fullRange[1] = newRange[1];
    flip->SetStencil(cylinder->GetOutput());
    std::vector<std::string> components;
    vtksys::SystemTools::SplitPath(fileName, components);
    cout << "extracting circular ROI from " << components.back() << endl;
    }
  flip->SetInputConnection(reader->GetOutputPort());
  flip->SetResliceAxesDirectionCosines(-1,0,0, 0,-1,0, 0,0,1);
  flip->SetOutputOrigin(origin);
  flip->SetOutputExtent(extent);
  flip->Update();

  image = flip->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  // generate the matrix, but modify to use DICOM coords
  static double xyFlipMatrix[16] =
    { -1, 0, 0, 0,  0, -1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
  // correct for the flip that was done earlier
  vtkMatrix4x4::Multiply4x4(*reader->GetDirectionCosines()->Element,
                            xyFlipMatrix, *matrix->Element);
  // do the left/right, up/down dicom-to-minc transformation
  vtkMatrix4x4::Multiply4x4(xyFlipMatrix, *matrix->Element, *matrix->Element);
  matrix->Modified();

  return true;
}
