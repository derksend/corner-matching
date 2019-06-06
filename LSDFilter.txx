/*
 * Copyright (C) 1999-2011 Insight Software Consortium
 * Copyright (C) 2005-2017 Centre National d'Etudes Spatiales (CNES)
 *
 * This file is part of Orfeo Toolbox
 *
 *     https://www.orfeo-toolbox.org/
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef LSDFilter_txx
#define LSDFilter_txx
#include "LSDFilter.h"

namespace otb
{

  template<class TImageType>
  LSDFilter<TImageType>
  ::LSDFilter() : m_layer(NULL, false)
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);
  }

  template<class TImageType>
  LSDFilter<TImageType>
  ::~LSDFilter()
  {
  }

  template<class TImageType>
  typename LSDFilter<TImageType>::OutputVectorDataPointerType
  LSDFilter<TImageType>
  ::ProcessTile()
  {
    std::cout << "1" << "\n";
    //OGRLayer m_layer(NULL, false);
    
    OGRSpatialReference oSRS(this->GetInput()->GetProjectionRef().c_str());
    std::vector<std::string> options;
    std::cout << "2" << "\n";
    
    typename ogr::DataSource::Pointer ogrDS = ogr::DataSource::New("test.shp", otb::ogr::DataSource::Modes::Overwrite);
    m_layer = ogrDS->CreateLayer("lines", &oSRS, wkbMultiLineString, options);
    std::cout << "3" << "\n";
    
    OGRFieldDefn fieldWidth("width",OFTReal);
    OGRFieldDefn fieldP("p",OFTReal);
    OGRFieldDefn fieldLogNFA("-log(NFA)",OFTReal);
    std::cout << "4" << "\n";
    m_layer.CreateField(fieldWidth);
    m_layer.CreateField(fieldP);
    m_layer.CreateField(fieldLogNFA);
    
    //Sums calculation per label
    unsigned int currentProgress = 0;
    unsigned long totalNbSegments = 0;
    std::cout << "5" << "\n";
    typedef otb::ExtractROI<typename TImageType::PixelType, typename TImageType::PixelType> ExtractImageFilterType;
    typename ExtractImageFilterType::Pointer extract = ExtractImageFilterType::New();
    extract->SetInput(this->GetInput());
    extract->SetExtractionRegion(this->GetInput()->GetBufferedRegion());
    extract->Update();
    std::cout << "6" << "\n";

    unsigned long sizeX=this->GetInput()->GetBufferedRegion().GetSize()[0];
    unsigned long sizeY=this->GetInput()->GetBufferedRegion().GetSize()[1];

    double * image = new double[(size_t)(sizeX*sizeY)];
    std::cout << sizeX << " " << sizeY << "\n";

    ImageIteratorType it(extract->GetOutput(),extract->GetOutput()->GetLargestPossibleRegion());

    unsigned long idx = 0;
    for(it.GoToBegin();!it.IsAtEnd() && idx < sizeX*sizeY;++it,++idx)
      {
	image[idx] = it.Get();
      }
    std::cout << "8" << "\n";
    std::cout << image[0] << std::endl;

    // int N = 10;
    
    // for (unsigned i = 0; i < N; i++) {
    //   double scale =i/10.0+0.01;
    //   for (unsigned i = 0; i < N; i++) {
    // 	double sigma_coef=0.01+i/10.0;	  
    // 	for (unsigned i = 0; i < N; i++) {
    // 	  double ang_th=0.01+i/10.0;
	    
    // 	  int n = 0;
    // 	  int regXNotUsed = 0;
    // 	  int regYNotUsed = 0;
    // 	  double log_eps=0;
    // 	  double density_th=0.7;
    // 	  double n_bins=1024;
    // 	  double quant=0;
    // 	  double * segs = LineSegmentDetection(&n, image, sizeX, sizeY,
    // 					       scale, sigma_coef, quant, 
    // 					       ang_th, log_eps, density_th, 
    // 					       n_bins,NULL,&regXNotUsed,&regYNotUsed);
    // 	  std::cout <<scale << " " << sigma_coef << " " << quant << " " << ang_th << " " <<  "n=" << n << "\n";
    // 	  if (n > 0){
    // 	    std::cout << "STOP" << "\n";
    // 	  }
    // 	}
    //   }
    // }

    int n = 0;
    int regXNotUsed = 0;
    int regYNotUsed = 0;
    double scale =.11;
    double sigma_coef=.11;
    double quant=0;
    double ang_th=0.91;
    double log_eps=0;
    double density_th=0.7;
    double n_bins=1024;
    double * segs = LineSegmentDetection(&n, image, sizeX, sizeY,
					 scale, sigma_coef, quant, 
					 ang_th, log_eps, density_th, 
					 n_bins,NULL,&regXNotUsed,&regYNotUsed);
    std::cout << "n=" << n << "\n";

    //std::cout << "n=" << n << "\n";

    m_totalNbSegments+=n;

    for(unsigned int idx = 0; idx < n;++idx)
      {
	OGRMultiLineString newMultiLineString;
	OGRLineString      newLineString;

	itk::ContinuousIndex<double,2> a,b;
	a[0] = segs[idx * 7];
	a[1] = segs[idx * 7 + 1];
	b[0] = segs[idx * 7 + 2];
	b[1] = segs[idx * 7 + 3];

	double width = segs[idx * 7 + 4];
	double p = segs[idx * 7 + 5];
	double lognfa = segs[idx * 7 + 6];

	typename TImageType::PointType pa, pb;

	extract->GetOutput()->TransformContinuousIndexToPhysicalPoint(a,pa);
	extract->GetOutput()->TransformContinuousIndexToPhysicalPoint(b,pb);

	// We assume isotropic spacing here
	width *= vcl_abs(extract->GetOutput()->GetSpacing()[0]);
           
	newLineString.addPoint(pa[0],pa[1],0);
	newLineString.addPoint(pb[0],pb[1],0);

	newMultiLineString.addGeometry(&newLineString);

	otb::ogr::Feature newFeature(m_layer.GetLayerDefn());

	newFeature.ogr().SetGeometry(&newMultiLineString);
	newFeature.ogr().SetField("Width",width);
	newFeature.ogr().SetField("-log(NFA)",lognfa);
	newFeature.ogr().SetField("p",p);
	m_layer.CreateFeature(newFeature);
      }
    std::cout << "10" << "\n";

    m_layer.ogr().CommitTransaction();
    std::cout << m_layer.GetFeatureCount(true) << "\n";

    ogr::Layer::feature_iter<ogr::Feature> fi = m_layer.begin();
    typedef typename VectorDataType::DataNodeType  DataNodeType;
    OGRIOHelper::Pointer helper = OGRIOHelper::New();

    typename OutputVectorDataType::Pointer outData = OutputVectorDataType::New();
    typename DataNodeType::Pointer root = outData->GetDataTree()->GetRoot()->Get();

    typename DataNodeType::Pointer document = DataNodeType::New();
    document->SetNodeType(otb::DOCUMENT);
    // // Adding the layer to the data tree
    outData->GetDataTree()->Add(document, root);
    typename DataNodeType::Pointer folder = DataNodeType::New();
    folder->SetNodeType(otb::FOLDER);
    outData->GetDataTree()->Add(folder, document);
    for (;fi!=m_layer.end();fi++) {
      DataNodeType::Pointer node = DataNodeType::New();
      helper->ConvertGeometryToLineNode(fi->GetGeometry(), node) ;
      outData->GetDataTree()->Add(node, folder);
    }
    
    // Release memory
    // delete [] image;
    // free( (void *) segs );

    // // OGRLayer l(m_layer);

    // std::cout << "12" << "\n";
    // 
    
    // std::cout << "13" << "\n";

    // //outData->GetDataTree()->SetRoot(rootNode);

    // typename DataNodeType::Pointer folder = DataNodeType::New();
    // folder->SetNodeType(otb::FOLDER);
    // outData->GetDataTree()->Add(folder, document);
    // std::cout << "14" << "\n";
    // outData->SetProjectionRef(this->GetInput()->GetProjectionRef());
    std::cout << outData->Size() << "\n";
    return outData;
    
  }


} // end namespace otb
#endif
