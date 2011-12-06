/*=========================================================================

Program:   EpilepsyViewerData
Module:    EpilepsyViewerData.h

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

#ifndef __EpilepsyViewerData_h
#define __EpilepsyViewerData_h

#include <vtkSmartPointer.h>

class vtkImageData;
class vtkPolyData;
class vtkMatrix4x4;
class vtkTimerLog;

//! Handle all data IO and data processing for the EpilepsyViewer
class EpilepsyViewerData
{
public:
  EpilepsyViewerData();
  ~EpilepsyViewerData();

  //! Load everything that is available in the data directory.
  //! This includes the images, the registration matrix, etc.
  bool LoadFromDataDirectory(const char *dirname);

  //! Register the MR head image to the CT head image.
  //! Once this is done, the MRHeadMatrix will be updated so
  //! that it puts the MR image into the same patient coordinate
  //! system as the CT image.
  bool RegisterMRHeadToCTHead();

  //! Extract the brain from the MR head image.
  bool ExtractMRBrain();

  //! Extract the electrodes from the CT head image.
  //! The electrodes will be in the CT data coordinate system,
  //! not the CT patient coordinate system, so when they are
  //! dislayed the CT matrix must be used to show them in the
  //! patient coordinate system.
  bool ExtractCTElectrodes();

  //! Get the CT head image.
  vtkImageData *GetCTHeadImage();

  //! Get the CT head data range.
  void GetCTHeadAutoRange(double r[2]) {
    r[0] = this->CTHeadAutoRange[0]; r[1] = this->CTHeadAutoRange[1]; }
  void GetCTHeadFullRange(double r[2]) {
    r[0] = this->CTHeadFullRange[0]; r[1] = this->CTHeadFullRange[1]; }

  //! Get the CT head image matrix (data coords to patient coords).
  vtkMatrix4x4 *GetCTHeadMatrix();

  //! Get the MR head image.
  vtkImageData *GetMRHeadImage();

  //! Get the MR head data range.
  void GetMRHeadAutoRange(double r[2]) {
    r[0] = this->MRHeadAutoRange[0]; r[1] = this->MRHeadAutoRange[1]; }
  void GetMRHeadFullRange(double r[2]) {
    r[0] = this->MRHeadFullRange[0]; r[1] = this->MRHeadFullRange[1]; }

  //! Get the MR brain image.
  vtkImageData *GetMRBrainImage();

  //! Get the MR brain surface
  vtkPolyData *GetMRBrainSurface();

  //! Get the MR head image matrix (data coords to patient coords).
  vtkMatrix4x4 *GetMRHeadMatrix();

  //! Get the electrodes.
  vtkPolyData *GetCTElectrodes();

private:
  //! Read a MINC image into the supplied data container.
  //! The direction cosines matrix will also be read from the file.
  //! The data will be converted to DICOM row and slice ordering,
  //! and the direction cosines will be modified so that the provide
  //! a conversion to the DICOM patient coordinate system.
  bool ReadMINCImage(
    vtkImageData *data, vtkMatrix4x4 *matrix,
    double fullRange[2], double autoRange[2],
    const char *filename);

  //! Read a DICOM series into the supplied data container.
  //! The DICOM data must not be compressed, or it cannot be read.
  //! The position and orientation will be converted into a 4x4 matrix.
  //! You must provide the name of a directory that contains only the
  //! desired image series and no other files.
  bool ReadDICOMSeries(
    vtkImageData *data, vtkMatrix4x4 *matrix,
    double fullRange[2], double autoRange[2],
    const char *directoryname);

  //! Report an error with printf formatting.
  void ReportError(const char *format, ...);

  vtkSmartPointer<vtkImageData> CTHeadImage;
  vtkSmartPointer<vtkImageData> MRHeadImage;
  vtkSmartPointer<vtkMatrix4x4> CTHeadMatrix;
  vtkSmartPointer<vtkMatrix4x4> MRHeadMatrix;
  double CTHeadAutoRange[2];
  double CTHeadFullRange[2];
  double MRHeadAutoRange[2];
  double MRHeadFullRange[2];

  vtkSmartPointer<vtkImageData> MRBrainImage;
  vtkSmartPointer<vtkPolyData> MRBrainSurface;

  vtkSmartPointer<vtkPolyData> CTElectrodes;

  vtkSmartPointer<vtkTimerLog> Timer;
};

#endif
