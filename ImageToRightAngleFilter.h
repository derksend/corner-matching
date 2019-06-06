#ifndef ImageToRightAngleFilter_h
#define ImageToRightAngleFilter_h

#include "LSDFilter.h"
#include "otbVectorDataToRightAngleVectorDataFilter.h"
#include "otbPersistentImageToVectorDataFilter.h"
#include "otbPersistentFilterStreamingDecorator.h"
#include "otbGridResampleImageFilter.h"
#include "otbConcatenateVectorDataFilter.h"
#include "otbImageFileWriter.h"
#include "otbExtractROI.h"
#include "otbOGRDataSourceWrapper.h"
#include "otbBandMathImageFilter.h"
// #include "ImageToRightAngleFilter.h"
#include "otbVectorDataTransformFilter.h"
#include "itkAffineTransform.h"

#include <sstream>
#include <unordered_set>
#include <string>
namespace otb
{

template <class TImageType>
class ImageToRightAngleFilter
  : public otb::PersistentImageToVectorDataFilter<TImageType, typename otb::LSDFilter<TImageType>::VectorDataType >
{
public:
  /** Standard Self typedef */
  typedef ImageToRightAngleFilter                                  Self;
  typedef typename otb::LSDFilter<TImageType>::VectorDataType VectorDataType;
  typedef PersistentImageToVectorDataFilter<TImageType, VectorDataType> Superclass;
  typedef itk::SmartPointer<Self>                                             Pointer;
  typedef itk::SmartPointer<const Self>                                       ConstPointer;

  typedef otb::LSDFilter<TImageType>     LSDType;
  typedef otb::VectorDataToRightAngleVectorDataFilter<VectorDataType> RightAngleFilterType;

  typedef typename Superclass::InputImageType              InputImageType;
  typedef typename Superclass::InputImagePointer           InputImagePointerType;

  typedef typename Superclass::OutputVectorDataType        OutputVectorDataType;
  typedef typename Superclass::OutputVectorDataPointerType OutputVectorDataPointerType;

  typedef typename OutputVectorDataType::DataTreeType DataTreeType;
  
  typedef itk::PreOrderTreeIterator<typename VectorDataType::DataTreeType> TreeIteratorType;
  typedef otb::ConcatenateVectorDataFilter<VectorDataType> ConcatenateVectorDataFilterType;
  typedef otb::BandMathImageFilter<TImageType> BandMathFilterType;
  typedef typename OutputVectorDataType::DataNodeType DataNodeType;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ImageToRightAngleFilter, PersistentImageToVectorDataFilter);

  itkGetMacro(DistanceThreshold, double);
  itkSetMacro(DistanceThreshold, double);

  itkGetMacro(AngleThreshold, double);
  itkSetMacro(AngleThreshold, double);

  itkGetMacro(MatchThreshold, double);
  itkSetMacro(MatchThreshold, double);

  itkGetMacro(NumberOfMatchingCorners, unsigned);
  itkGetMacro(NumberOfCornersInInput, unsigned);
  itkGetMacro(NumberOfCornersInInput2, unsigned);

  itkSetMacro(Tempdir,std::string);

  itkSetMacro(Scale, double);
  itkSetMacro(Sigmacoef, double);
  itkSetMacro(Quant, double);
  itkSetMacro(Angth, double);
  itkSetMacro(Logeps, double);
  itkSetMacro(Densityth, double);
  
  void SetInput2(TImageType* inputImage);
  const TImageType* GetInput2();
  
protected:
  ImageToRightAngleFilter();

  ~ImageToRightAngleFilter() ITK_OVERRIDE;

  OutputVectorDataPointerType ExtractCorners(const TImageType* inputImage);
  OutputVectorDataPointerType SimplifyGeometry(OutputVectorDataPointerType corners);
  void FixNodeIds(OutputVectorDataPointerType corners);

  //void GenerateInputRequestedRegion() ITK_OVERRIDE;

private:
  ImageToRightAngleFilter(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented

  OutputVectorDataPointerType ProcessTile() ITK_OVERRIDE;

  double m_DistanceThreshold;
  double m_AngleThreshold;
  double m_MatchThreshold;
  unsigned m_NumberOfMatchingCorners;
  unsigned m_NumberOfCornersInInput;
  unsigned m_NumberOfCornersInInput2;
  std::string m_Tempdir;
  double m_Scale;
  double m_Sigmacoef;
  double m_Quant;
  double m_Angth;
  double m_Logeps;
  double m_Densityth;
};

template <class TImageType>
class StreamingImageToRightAngleFilter
{
public:

  // define the ImageToRightAngleFilter template
  typedef ImageToRightAngleFilter<TImageType>
    ImageToRightAngleFilterType;

  typedef typename ImageToRightAngleFilterType::InputImageType
      InputImageType;
  typedef typename ImageToRightAngleFilterType::OutputVectorDataType
      OutputVectorDataType;

  // typedef for streaming capable filter
  typedef PersistentFilterStreamingDecorator<ImageToRightAngleFilterType>
    FilterType;

};

}

#ifndef OTB_MANUAL_INSTANTIATION
#include "ImageToRightAngleFilter.txx"
#endif

#endif



 
