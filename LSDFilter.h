#ifndef LSDFilter_h
#define LSDFilter_h

extern "C"
{
#include "lsd.h"
}
#include "otbVectorDataToRightAngleVectorDataFilter.h"
#include "otbPersistentImageToVectorDataFilter.h"
#include "otbPersistentFilterStreamingDecorator.h"
#include "otbGridResampleImageFilter.h"
#include "otbConcatenateVectorDataFilter.h"
#include "otbImageFileWriter.h"
#include "otbExtractROI.h"
//#include "otbVectorData.h"
// #include "otbOGRLayerWrapper.h"
// #include "otbOGRFeatureWrapper.h"
// #include "otbOGRFieldWrapper.h"
#include "otbOGRDataSourceWrapper.h"
#include "otbBandMathImageFilter.h"
#include "LSDFilter.h"
#include "otbVectorDataTransformFilter.h"
#include "itkAffineTransform.h"
#include "otbImage.h"
#include "otbOGRIOHelper.h"

#include <sstream>

namespace otb
{

template <class TImageType>
class LSDFilter
  : public otb::PersistentImageToVectorDataFilter<TImageType, otb::VectorData<double> >
{
public:
  /** Standard Self typedef */
  typedef LSDFilter                                  Self;
  typedef typename otb::VectorData<double> VectorDataType;
  typedef PersistentImageToVectorDataFilter<TImageType, VectorDataType> Superclass;
  typedef itk::SmartPointer<Self>                                             Pointer;
  typedef itk::SmartPointer<const Self>                                       ConstPointer;

  typedef typename Superclass::InputImageType              InputImageType;
  typedef typename Superclass::InputImagePointer           InputImagePointerType;

  typedef typename Superclass::OutputVectorDataType        OutputVectorDataType;
  typedef typename Superclass::OutputVectorDataPointerType OutputVectorDataPointerType;

  typedef itk::ImageRegionConstIterator<InputImageType> ImageIteratorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(LSDFilter, PersistentImageToVectorDataFilter);

  itkGetMacro(totalNbSegments, unsigned);
  itkGetMacro(layer, otb::ogr::Layer);
  
  
protected:
  LSDFilter();

  ~LSDFilter() ITK_OVERRIDE;

  OutputVectorDataPointerType ExtractCorners(const TImageType* inputImage);

  //void GenerateInputRequestedRegion() ITK_OVERRIDE;

private:
  LSDFilter(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented

  OutputVectorDataPointerType ProcessTile() ITK_OVERRIDE;

  unsigned m_totalNbSegments;
  otb::ogr::Layer m_layer;
  
};

}

#ifndef OTB_MANUAL_INSTANTIATION
#include "LSDFilter.txx"
#endif

#endif



 
