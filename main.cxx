#include "otbImageFileReader.h"
#include "otbStreamingLineSegmentDetector.h"
#include "otbAlphaBlendingFunctor.h"
#include "itkBinaryFunctorImageFilter.h"
#include "otbImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "otbVectorDataFileWriter.h"
#include "otbVectorDataFileReader.h"
#include "otbVectorDataToRightAngleVectorDataFilter.h"
#include "ImageToRightAngleFilter.h"
#include "otbVectorImage.h"
#include <unordered_set>
int main(int argc, char * argv[])
{

  if (argc < 14) {
    std::cout << "Usage : input1 input2 output angleThreshold distanceThreshold nbX matchThreshold tmpdir scale sigmacoef quant angth logeps densityth\n";
    return 0;
  }
  const char * input1             = argv[1];
  const char * input2             = argv[2];
  const char * outfname           = argv[3];
  double       angleThreshold     = atof(argv[4]);
  double       distanceThreshold  = atof(argv[5]);
  unsigned     nbX                = atoi(argv[6]);
  double       matchThreshold     = atof(argv[7]);
  std::string  tmpdir(argv[8]);
  double       scale              = atof(argv[9]);
  double       sigmacoef          = atof(argv[10]);
  double       quant              = atof(argv[11]);
  double       angth              = atof(argv[12]);
  double       logeps             = atof(argv[13]);
  double       densityth          = atof(argv[14]);

  GDALSetCacheMax(0);
  typedef double                                                         PrecisionType;
  typedef unsigned                                                       ComponentType;
  
  typedef otb::Image<ComponentType>                                      ImageType;
  typedef otb::ImageFileReader<ImageType>                                ImageReaderType;
    
  typedef otb::VectorData<PrecisionType>                                 VectorDataType;
  typedef otb::VectorDataFileWriter<VectorDataType>                      WriterType;

  ImageReaderType::Pointer reader1 = ImageReaderType::New();
  reader1->SetFileName(input1);
  reader1->GenerateOutputInformation();
  // reader1->Update();

  ImageReaderType::Pointer reader2 = ImageReaderType::New();
  reader2->SetFileName(input2);
  reader2->GenerateOutputInformation();
  // reader2->Update();

   
  
  typedef typename otb::StreamingImageToRightAngleFilter<ImageType>::FilterType ImageToRightAngleFilterType;
  
  ImageToRightAngleFilterType::Pointer rightAngleFilter  = ImageToRightAngleFilterType::New();

  rightAngleFilter->GetFilter()->SetInput(reader1->GetOutput());
  rightAngleFilter->GetFilter()->SetInput2(reader2->GetOutput());

  rightAngleFilter->GetStreamer()->SetNumberOfDivisionsTiledStreaming(nbX*nbX);

  rightAngleFilter->GetFilter()->SetAngleThreshold(angleThreshold);
  rightAngleFilter->GetFilter()->SetDistanceThreshold(distanceThreshold);
  rightAngleFilter->GetFilter()->SetMatchThreshold(matchThreshold);
  rightAngleFilter->GetFilter()->SetTempdir(tmpdir);
  rightAngleFilter->GetFilter()->SetScale(scale);
  rightAngleFilter->GetFilter()->SetSigmacoef(sigmacoef);
  rightAngleFilter->GetFilter()->SetQuant(quant);
  rightAngleFilter->GetFilter()->SetAngth(angth);
  rightAngleFilter->GetFilter()->SetLogeps(logeps);
  rightAngleFilter->GetFilter()->SetDensityth(densityth);

  rightAngleFilter->Update();
  // std::cout << "Update finished" << "\n";
  // std::cout << "After" << rightAngleFilter->GetFilter()->GetOutputVectorData()->Size() << "\n";

  WriterType::Pointer writer = WriterType::New();
  
  writer->SetInput(rightAngleFilter->GetFilter()->GetOutputVectorData());
  
  writer->SetFileName(outfname);

  writer->Update();
  //std::cout << rightAngleFilter->GetFilter()->GetNumberOfMatchingCorners() << "\n";
  //std::cout << rightAngleFilter->GetFilter()->GetNumberOfCornersInInput() << "\n";

  std::cout << rightAngleFilter->GetFilter()->GetNumberOfCornersInInput() << " " << rightAngleFilter->GetFilter()->GetNumberOfCornersInInput2() << " " << (float)rightAngleFilter->GetFilter()->GetNumberOfMatchingCorners()/rightAngleFilter->GetFilter()->GetNumberOfCornersInInput() << "\n";
  
  return EXIT_SUCCESS;
}
