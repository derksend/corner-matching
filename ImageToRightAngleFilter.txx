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

#ifndef ImageToRightAngleFilter_txx
#define ImageToRightAngleFilter_txx
#include "ImageToRightAngleFilter.h"

namespace otb
{

  template<class TImageType>
  ImageToRightAngleFilter<TImageType>
  ::ImageToRightAngleFilter() : m_DistanceThreshold(500),
				m_AngleThreshold(0.8),
				m_MatchThreshold(500),
				m_NumberOfMatchingCorners(0),
				m_NumberOfCornersInInput(0),
				m_NumberOfCornersInInput2(0),
				m_Scale(0.8),
				m_Sigmacoef(0.6),
				m_Quant(2),
				m_Angth(22.5),
				m_Logeps(0),
				m_Densityth(0.7)
  {
    this->SetNumberOfRequiredInputs(2);
  }

  template<class TImageType>
  ImageToRightAngleFilter<TImageType>
  ::~ImageToRightAngleFilter()
  {
  }
  template<class TImageType>
  void
  ImageToRightAngleFilter<TImageType>
  ::SetInput2(TImageType* inputImage)
  {
    this->SetNthInput(1,inputImage);
  }

  template<class TImageType>
  const TImageType*
  ImageToRightAngleFilter<TImageType>
  ::GetInput2()
  {
    return this->GetInput(1);
  }

  
  template<class TImageType>
  typename ImageToRightAngleFilter<TImageType>::OutputVectorDataPointerType
  ImageToRightAngleFilter<TImageType>
  ::ExtractCorners(const TImageType* inputImage)
  {
    typename ConcatenateVectorDataFilterType::Pointer conc = ConcatenateVectorDataFilterType::New();
    itk::ImageRegionConstIterator<TImageType> it(inputImage,inputImage->GetLargestPossibleRegion());
    std::cout << "start" << "\n";
    std::unordered_set<unsigned> labelList;
    for(it.GoToBegin();!it.IsAtEnd() ;++it)
      {
	labelList.insert(it.Get());
      }
    std::vector<unsigned> m_ClassList(labelList.begin(),labelList.end());
    std::sort(m_ClassList.begin(),m_ClassList.end());
    for (auto a : m_ClassList) {
      std::cout << a << "\n";
    }
    if (m_ClassList.size() > 1){
    for (std::vector<unsigned>::const_iterator i = m_ClassList.begin(); i != m_ClassList.end(); i++) {
      typename BandMathFilterType::Pointer bm = BandMathFilterType::New();
      bm->SetNthInput(0, const_cast<TImageType*>(inputImage));
      std::stringstream s;
      s << "b1==" << *i << "?255:1";
      bm->SetExpression(s.str());
      bm->Update();

    
      OGRSpatialReference oSRS(this->GetInput()->GetProjectionRef().c_str());
      std::vector<std::string> options;
      std::stringstream dataSetName;
      dataSetName << m_Tempdir << "lines_class"<< *i << ".shp";
      typename ogr::DataSource::Pointer ogrDS = ogr::DataSource::New(dataSetName.str(), otb::ogr::DataSource::Modes::Overwrite);
      ogr::Layer m_layer = ogrDS->CreateLayer("lines", &oSRS, wkbMultiLineString, options);
    
      OGRFieldDefn fieldWidth("width",OFTReal);
      OGRFieldDefn fieldP("p",OFTReal);
      OGRFieldDefn fieldLogNFA("-log(NFA)",OFTReal);
      m_layer.CreateField(fieldWidth);
      m_layer.CreateField(fieldP);
      m_layer.CreateField(fieldLogNFA);
    
      //Sums calculation per label
      unsigned int currentProgress = 0;
      unsigned long totalNbSegments = 0;

      unsigned long sizeX=bm->GetOutput()->GetBufferedRegion().GetSize()[0];
      unsigned long sizeY=bm->GetOutput()->GetBufferedRegion().GetSize()[1];

      double * image = new double[(size_t)(sizeX*sizeY)];
      typedef itk::ImageRegionConstIterator<TImageType> ImageIteratorType;

      ImageIteratorType it(bm->GetOutput(),bm->GetOutput()->GetBufferedRegion());

      unsigned long idx = 0;
      for(it.GoToBegin();!it.IsAtEnd() && idx < sizeX*sizeY;++it,++idx)
	{
	  image[idx] = it.Get();
	}


      int n = 0;
      int regXNotUsed = 0;
      int regYNotUsed = 0;
      double n_bins=1024;
      double * segs = LineSegmentDetection(&n, image, sizeX, sizeY,
					   m_Scale, m_Sigmacoef, m_Quant, 
					   m_Angth, m_Logeps, m_Densityth, 
					   n_bins,NULL,&regXNotUsed,&regYNotUsed);

      totalNbSegments+=n;

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

	  bm->GetOutput()->TransformContinuousIndexToPhysicalPoint(a,pa);
	  bm->GetOutput()->TransformContinuousIndexToPhysicalPoint(b,pb);

	  // We assume isotropic spacing here
	  width *= vcl_abs(bm->GetOutput()->GetSpacing()[0]);
           
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

      m_layer.ogr().CommitTransaction();

      // Release memory
      delete [] image;
      free( (void *) segs );

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
	typename DataNodeType::Pointer node = DataNodeType::New();
	helper->ConvertGeometryToLineNode(fi->GetGeometry(), node) ;
	outData->GetDataTree()->Add(node, folder);
      }
      conc->AddInput(outData);
    }
    conc->Update();
    typename RightAngleFilterType::Pointer rightAngleFilter  = RightAngleFilterType::New();
    rightAngleFilter->SetInput(conc->GetOutput());
    rightAngleFilter->SetAngleThreshold(m_AngleThreshold);
    rightAngleFilter->SetDistanceThreshold(m_DistanceThreshold);
    rightAngleFilter->Update();
    return rightAngleFilter->GetOutput();
    }
    else{
      return VectorDataType::New();
    }
  }

  template<class TImageType>
  void
  ImageToRightAngleFilter<TImageType>
  ::FixNodeIds(OutputVectorDataPointerType corners){
    TreeIteratorType i1(corners->GetDataTree());
    i1.GoToBegin();
    unsigned i = 0;
    for (; !i1.IsAtEnd(); i1++) {
      if (i1.Get()->IsPointFeature()){
	std::stringstream ss;
	ss << i1.Get()->GetNodeId() << i++;
	i1.Get()->SetNodeId(ss.str());
      }
    }
  }
  
  template<class TImageType>
  typename ImageToRightAngleFilter<TImageType>::OutputVectorDataPointerType
  ImageToRightAngleFilter<TImageType>
  ::SimplifyGeometry(OutputVectorDataPointerType corners)
  {
    std::cout << "start simplification" << "\n";
    std::cout << "Initial corners : " << corners->Size()-3 << "\n";

    //Make node ID unique for faster computation later
    FixNodeIds(corners);

    typename OutputVectorDataType::Pointer outData = OutputVectorDataType::New();
    // Retrieving root node
    typename DataNodeType::Pointer root = outData->GetDataTree()->GetRoot()->Get();
    // Create the document node
    typename DataNodeType::Pointer document = DataNodeType::New();
    document->SetNodeType(otb::DOCUMENT);
    // Adding the layer to the data tree
    outData->GetDataTree()->Add(document, root);
    // Create the folder node
    typename DataNodeType::Pointer folder = DataNodeType::New();
    folder->SetNodeType(otb::FOLDER);
    outData->GetDataTree()->Add(folder, document);
    outData->SetProjectionRef(this->GetInput()->GetProjectionRef());
    
    // Itterate on the vector data
    TreeIteratorType i1(corners->GetDataTree());      // Reference
    i1.GoToBegin();
    TreeIteratorType i2(corners->GetDataTree());      // Reference
    std::unordered_set<std::string> usedPoints;

    for (; !i1.IsAtEnd(); i1++) {


      if (!i1.Get()->IsPointFeature() || usedPoints.find(i1.Get()->GetNodeId()) != usedPoints.end())
    	{
    	  continue; // do not process if it's not a point
    	}
      i2.GoToBegin();
      unsigned numberOfNeighbors=0;
      std::vector<typename DataNodeType::PointType> neighbors;
      for (; !i2.IsAtEnd(); i2++) {
    	if (!i2.Get()->IsPointFeature() || usedPoints.find(i2.Get()->GetNodeId()) != usedPoints.end())
    	  {
    	    continue; // do not process if it's not a point
    	  }

    	if (i1.Get()->GetPoint().EuclideanDistanceTo(i2.Get()->GetPoint()) < m_MatchThreshold ){
	  numberOfNeighbors++;
	  neighbors.push_back(i2.Get()->GetPoint());
	  usedPoints.insert(i2.Get()->GetNodeId());
	}
      }
      if (numberOfNeighbors > 1) {

	double p[2];
	p[0]=0;
	p[1]=0;
	for (unsigned i = 0; i < numberOfNeighbors; i++) {
	  p[0]+=neighbors[i][0]/numberOfNeighbors;
	  p[1]+=neighbors[i][1]/numberOfNeighbors;
	}
	typename DataNodeType::PointType centroid(p);

	usedPoints.insert(i1.Get()->GetNodeId());
	typename DataNodeType::Pointer CurrentGeometry = DataNodeType::New();
	CurrentGeometry->SetNodeId("FEATURE_POINT");
	CurrentGeometry->SetNodeType(otb::FEATURE_POINT);
	CurrentGeometry->SetPoint(centroid);
	outData->GetDataTree()->Add(CurrentGeometry, folder);
      }
      else{
	usedPoints.insert(i1.Get()->GetNodeId());
	typename DataNodeType::Pointer CurrentGeometry = DataNodeType::New();
	CurrentGeometry->SetNodeId("FEATURE_POINT");
	CurrentGeometry->SetNodeType(otb::FEATURE_POINT);
	CurrentGeometry->SetPoint(i1.Get()->GetPoint());
	outData->GetDataTree()->Add(CurrentGeometry, folder);
      }
    }
    return outData;
  }
  
  template<class TImageType>
  typename ImageToRightAngleFilter<TImageType>::OutputVectorDataPointerType
  ImageToRightAngleFilter<TImageType>
  ::ProcessTile()
  {
    //std::cout << "Start new tile" << "\n";
    // Apply an ExtractImageFilter to avoid problems with filters asking for the LargestPossibleRegion
    typedef otb::ExtractROI<typename TImageType::PixelType, typename TImageType::PixelType> ExtractImageFilterType;
    typename ExtractImageFilterType::Pointer extract = ExtractImageFilterType::New();
    extract->SetInput(this->GetInput());
    extract->SetExtractionRegion(this->GetInput()->GetBufferedRegion());
    extract->Update();

    OutputVectorDataPointerType corners1 = ExtractCorners(extract->GetOutput());
    corners1 = SimplifyGeometry(corners1);
    std::cout << "After simplification: " << corners1->Size()-3 << "\n";
    m_NumberOfCornersInInput+=corners1->Size()-3; //-3 for data root nodes
	
    //Extract corners from image 2 with buffered region from input 1
    extract->SetInput(this->GetInput2());
    extract->SetExtractionRegion(this->GetInput()->GetBufferedRegion());
    extract->Update();
	
    OutputVectorDataPointerType corners2 = ExtractCorners(extract->GetOutput());

    corners2 = SimplifyGeometry(corners2);
    m_NumberOfCornersInInput2+=corners2->Size()-3;
    std::cout << "After simplification: " << corners2->Size()-3 << "\n";

    typedef typename VectorDataType::DataNodeType  DataNodeType;
    typedef typename VectorDataType::PointType     PointType;

    typename OutputVectorDataType::Pointer outData = OutputVectorDataType::New();
    // Retrieving root node
    typename DataNodeType::Pointer root = outData->GetDataTree()->GetRoot()->Get();
    // Create the document node
    typename DataNodeType::Pointer document = DataNodeType::New();
    document->SetNodeType(otb::DOCUMENT);
    // Adding the layer to the data tree
    outData->GetDataTree()->Add(document, root);
    // Create the folder node
    typename DataNodeType::Pointer folder = DataNodeType::New();
    folder->SetNodeType(otb::FOLDER);
    outData->GetDataTree()->Add(folder, document);
    outData->SetProjectionRef(this->GetInput()->GetProjectionRef());
    
    // Itterate on the vector data
    TreeIteratorType i1(corners1->GetDataTree());      // Reference
    i1.GoToBegin();
    TreeIteratorType i2(corners2->GetDataTree());      // Reference

    for (; !i1.IsAtEnd(); i1++) {
      if (!i1.Get()->IsPointFeature())
	{
	  continue; // do not process if it's not a point
	}
      i2.GoToBegin();
      for (; !i2.IsAtEnd(); i2++) {
	if (!i2.Get()->IsPointFeature())
	  {
	    continue; // do not process if it's not a point
	  }
	if (i1.Get()->GetPoint().EuclideanDistanceTo(i2.Get()->GetPoint()) < m_MatchThreshold){
	  m_NumberOfMatchingCorners++;
	  typename DataNodeType::Pointer CurrentGeometry = DataNodeType::New();
	  CurrentGeometry->SetNodeId("FEATURE_POINT");
	  CurrentGeometry->SetNodeType(otb::FEATURE_POINT);
	  CurrentGeometry->SetPoint(i1.Get()->GetPoint());
	  outData->GetDataTree()->Add(CurrentGeometry, folder);
	  break;
	}
      }
    }
    return outData;
  }


} // end namespace otb
#endif
