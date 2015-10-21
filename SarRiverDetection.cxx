/*========================================================================= */

// CRIM : Centre de recherche en informatique de Montréal

// Équipe Télédetection pour les catastrophes majeures (TCM).

// Programme : Extraction des rivières à partir des images Radar.

// Auteur : Moslem Ouled Sghaier

// Version : 0

/*========================================================================= */

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "otbSFSTexturesImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include <itkGrayscaleMorphologicalOpeningImageFilter.h>
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFlatStructuringElement.h"

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iostream>
#include "itkVector.h"
#include "itkBresenhamLine.h"
#include <vector>
#include <math.h>
#include <ctime>


#include <cstdlib>

// Seuillage
#include "itkBinaryThresholdImageFilter.h"

//Bresenham
#include "itkBresenhamLine.h"

// Image  to LabelMap and LabelImage
#include "itkImageRegionIterator.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include <itkLabelMapToBinaryImageFilter.h>
#include "otbWrapperApplication.h" 
#include "otbWrapperApplicationRegistry.h"
#include "otbWrapperApplicationFactory.h"
#include "otbWrapperTags.h"

#include "itkConnectedComponentImageFilter.h"
#include "otbPersistentVectorizationImageFilter.h"

//Utils
#include "itksys/SystemTools.hxx"
#include "itkListSample.h"

// Elevation handler
#include "otbWrapperElevationParametersHandler.h"

// Image thining
#include "itkBinaryThinningImageFilter.h"


// ceci sera utile pour les données vectorielles
#include "otbVectorData.h" 
#include "otbVectorDataFileReader.h" 
#include "otbVectorDataFileWriter.h"

#include "itkPreOrderTreeIterator.h" 
#include "otbObjectList.h" 
#include "otbPolygon.h"


//using namespace std;

extern "C" {
	#include "pde_toolbox_bimage.h"
	#include "pde_toolbox_defs.h"
//	#include "pde_toolbox_LSTB.h"
//	#include "ImageMagickIO.h"
}

#include "pathopenclose.h"


namespace otb
{
namespace Wrapper
{

class SarRiverDetection : public Application
{

public:
	typedef SarRiverDetection Self;
    typedef Application              Superclass;
    typedef itk::SmartPointer<Self>       Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

#define PI 3.14159265


/** Standard macro */
    itkNewMacro(Self);
    itkTypeMacro(SarRiverDetection, otb::Application);

typedef unsigned char CharPixelType; // IO
typedef otb::Image<CharPixelType, 2> CharImageType;


private:

void DoInit()
{

SetName("SarRiverDetection"); // Nécessaire
SetDocName("SarRiverDetection");
SetDocLongDescription("Un simple module pour l'extraction des rivières à partir des images radar");
SetDocLimitations("Les autres paramètres seront ajoutés plus tard");
SetDocAuthors("Moslem Ouled Sghaier");


AddParameter(ParameterType_InputImage,"in", "Input Image");
SetParameterDescription("in", "The input image");
AddParameter(ParameterType_OutputVectorData,"out", "Output Image");
SetParameterDescription("out","The output image");

}

void DoUpdateParameters()
{
	// Nothing to do here : all parameters are independent
}


void DoExecute()
{
  

  clock_t start, stop;
  start = clock();

   /////////////////////////////////////////////////////////////////////// déclaration des structures globales

typedef unsigned char CharPixelType; // IO
const unsigned int Dimension = 2;
typedef otb::Image<CharPixelType, Dimension> CharImageType;

typedef itk::RescaleIntensityImageFilter<UInt32ImageType,CharImageType> RescaleFilter0;
RescaleFilter0::Pointer rescale0 = RescaleFilter0::New();
rescale0->SetInput(GetParameterUInt32Image("in"));
rescale0->UpdateLargestPossibleRegion();


  // Application du SFS-SD

  typedef otb::SFSTexturesImageFilter<CharImageType, CharImageType> SFSFilterType;
  SFSFilterType::Pointer filter   = SFSFilterType::New();
  filter->SetSpectralThreshold(15);
  filter->SetSpatialThreshold(100); 
  filter->SetNumberOfDirections(20); 
  filter->SetRatioMaxConsiderationNumber(5); 
  filter->SetAlpha(1.00);

  filter->SetInput(rescale0->GetOutput());
  filter->Update();

  // Filtrer les bords de l'image

  unsigned int x= filter->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
  unsigned int y= filter->GetOutput()->GetLargestPossibleRegion().GetSize()[1];

  for (unsigned int i=0 ; i< 2;i++ )
  for (unsigned int j=0 ; j< y;j++)
  {{
  CharImageType::IndexType pixelIndex; 
  pixelIndex[0] = i;   // x position 
  pixelIndex[1] = j;   // y position
  filter->GetOutput()->SetPixel(pixelIndex,0);}}

  for (unsigned int i=0 ; i< x;i++ )
  for (unsigned int j=0 ; j< 2;j++)
  {{
  CharImageType::IndexType pixelIndex; 
  pixelIndex[0] = i;   // x position 
  pixelIndex[1] = j;   // y position
  filter->GetOutput()->SetPixel(pixelIndex,0);}}

  for (unsigned int i=x-2 ; i< x;i++ )
  for (unsigned int j=0 ; j< y;j++)
  {{
  CharImageType::IndexType pixelIndex; 
  pixelIndex[0] = i;   // x position 
  pixelIndex[1] = j;   // y position
  filter->GetOutput()->SetPixel(pixelIndex,0);}}

  for (unsigned int i=0 ; i< x;i++ )
  for (unsigned int j=y-2 ; j< y;j++)
  {{
  CharImageType::IndexType pixelIndex; 
  pixelIndex[0] = i;   // x position 
  pixelIndex[1] = j;   // y position
  filter->GetOutput()->SetPixel(pixelIndex,0);}}

  // Application de l'ouverture morphologique (opening)

  typedef itk::FlatStructuringElement<2>  StructuringElementType;
  typedef itk::GrayscaleMorphologicalOpeningImageFilter<CharImageType,CharImageType,StructuringElementType>  GrayscaleDilateImageFilterType;
  GrayscaleDilateImageFilterType::Pointer openingFilter = GrayscaleDilateImageFilterType::New();

  openingFilter->SetInput(filter->GetOutput());
  StructuringElementType::RadiusType elementRadius;
  elementRadius.Fill(3);
  StructuringElementType structuringElement = StructuringElementType::Box(elementRadius);
  openingFilter->SetKernel(structuringElement);
  openingFilter->GetSafeBorder();
  openingFilter->Update();

  // Fin de l'ouverture

  // Application de l'ouverture morhologique en utilisant un élément structurant sous la forme d'un chemin

  int   i, L, K;

  L=100; K=1;

  unsigned int nx = openingFilter->GetOutput()->GetLargestPossibleRegion().GetSize()[0];//input_bimage->dim->buf[0];
  unsigned int ny = openingFilter->GetOutput()->GetLargestPossibleRegion().GetSize()[1];;//input_bimage->dim->buf[1];
  unsigned int num_pixels = nx*ny;

  std::cout << num_pixels << std::endl;

  PATHOPEN_PIX_TYPE * input_image = new PATHOPEN_PIX_TYPE[nx * ny];
  PATHOPEN_PIX_TYPE * output_image = new PATHOPEN_PIX_TYPE[nx * ny];
  
  // Convert intermediate float to PATHOPEN_PIX_TYPE (unsigned char)
	for (i = 0; i < num_pixels; ++i) {
		CharImageType::IndexType index = {i % nx, i / nx};
		input_image[i] = openingFilter->GetOutput()->GetPixel(index);//static_cast<PATHOPEN_PIX_TYPE>(input_bimage->buf[i]);
	}

	std::cout << "Calling pathopen()" << std::endl;
        
	pathopen(
            input_image, // The input image //
            nx, ny,	 // Image dimensions //
            L,		 // The threshold line length //
            K,		 // The maximum number of gaps in the path //
            output_image // Output image //
            ); 

	for (i = 0; i < num_pixels; ++i) {
		CharImageType::IndexType index = {i % nx, i / nx};
		openingFilter->GetOutput()->SetPixel(index,  static_cast<CharImageType::PixelType>(output_image[i]));
	}

  // Fin de la morphologie mathématique

  // Seuillage pour afficher les zones homogènes

  typedef itk::BinaryThresholdImageFilter<CharImageType, CharImageType>  FilterType1;
  FilterType1::Pointer filter1 = FilterType1::New();
  filter1->SetInput(openingFilter->GetOutput());
  filter1->SetLowerThreshold(15); 
  filter1->Update();

  typedef itk::BinaryThresholdImageFilter<CharImageType, CharImageType>  FilterType2;
  FilterType2::Pointer filter2 = FilterType2::New();
  filter2->SetInput(openingFilter->GetOutput());
  filter2->SetLowerThreshold(15); 
  filter2->Update();
  
  // Fin du seuillage 
 
  // Faire le suivi des droites
  
  typedef itk::BinaryThinningImageFilter <CharImageType,CharImageType> BinaryThinningImageFilterType;
  BinaryThinningImageFilterType::Pointer binaryThinningImageFilter = BinaryThinningImageFilterType::New();
  binaryThinningImageFilter->SetInput(filter1->GetOutput());
  binaryThinningImageFilter->Update();

  // Rescale the image so that it can be seen (the output is 0 and 1, we want 0 and 255)
  typedef itk::RescaleIntensityImageFilter<CharImageType,CharImageType > RescaleType;
  RescaleType::Pointer rescaler = RescaleType::New();
  rescaler->SetInput(binaryThinningImageFilter->GetOutput());
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescaler->Update();
 
  std::vector<std::vector<CharImageType::IndexType>> V;
  V = Recherche(rescaler);

  for (int u=0;u< V.size();u++)
  {
	  for (int p=0;p<V.at(u).size();p++)
	  {
		  CharImageType::IndexType pixelIndex; 
          pixelIndex[0] = V.at(u).at(p)[0];   // x position 
          pixelIndex[1] = V.at(u).at(p)[1];   // y position
          filter1->GetOutput()->SetPixel(pixelIndex,255);
	  }
  } 
  
  // Fin du suivi
 
  // Filtrage par forme
  
  std::vector<int> v;
  
  typedef itk::BinaryImageToShapeLabelMapFilter<CharImageType> BinaryImageToShapeLabelMapFilterType;
  BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
  binaryImageToShapeLabelMapFilter->SetInput(filter1->GetOutput());
  binaryImageToShapeLabelMapFilter->FullyConnectedOn();
  binaryImageToShapeLabelMapFilter->Update();
  // The output of this filter is an itk::ShapeLabelMap, which contains itk::ShapeLabelObject's
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects 0." << std::endl;
  // Loop over all of the blobs
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++)
    {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    // Output the bounding box (an example of one possible property) of the ith region
	//std::cout << "Object " << i << " has bounding box " << labelObject->GetBoundingBox() << std::endl;
	if ( ( (labelObject->GetBoundingBox().GetSize()[0] * labelObject->GetBoundingBox().GetSize()[1]) / labelObject->GetNumberOfPixels() ) < 8 )
	{v.push_back(i);}
	}
  
  typedef itk::LabelObject<CharPixelType,2> LabelObjectType; 
  typedef itk::LabelMap<LabelObjectType>   LabelMapType;
  typedef itk::LabelMapToBinaryImageFilter<LabelMapType,CharImageType> L2BIType;
  typedef itk::BinaryImageToLabelMapFilter<CharImageType, LabelMapType> I2LType;

  I2LType::Pointer i2l = I2LType::New();
  i2l->SetInput(filter1->GetOutput());
  i2l->SetInputForegroundValue(255); 
  i2l->SetOutputBackgroundValue(0);
  i2l->FullyConnectedOn();
  i2l->Update();
  
  // The output of this filter is an itk::ShapeLabelMap, which contains itk::ShapeLabelObject's
  std::cout << "There are " << i2l->GetOutput()->GetNumberOfLabelObjects() << " objects." << std::endl;
  
  // Loop over all of the blobs
  int element=0;
  for(unsigned int i = 0; i < v.size(); i++)
  {
  I2LType::OutputImageType::LabelObjectType* labelObject = i2l->GetOutput()->GetNthLabelObject(v.at(i)+element);
  i2l->GetOutput()->RemoveLabelObject(labelObject); 
  element--;
  i2l->Update();
  }
  
  // Eliminer les segments de plus

  for (int u=0;u< V.size();u++)
  {
	  for (int p=1;p<V.at(u).size()-1;p++)
	  {
		  CharImageType::IndexType pixelIndex; 
          pixelIndex[0] = V.at(u).at(p)[0];   // x position 
          pixelIndex[1] = V.at(u).at(p)[1];   // y position
		  CharImageType::PixelType pixelValue = filter2->GetOutput()->GetPixel(pixelIndex);
		  if ( pixelValue == 0)
		  {i2l->GetOutput()->SetPixel(pixelIndex,0);}
	  }
  }
  i2l->Update();
   
  L2BIType::Pointer l2bi = L2BIType::New(); 
  l2bi->SetInput(i2l->GetOutput());
  l2bi->Update();
  
  // Fin du filtrage

  stop = clock();
  std::cout << "CPU time elapsed:" << ((double)stop-start)/CLOCKS_PER_SEC << std::endl;

  //Label the objects in a binary image

  typedef otb::Polygon<double>             PolygonType;
  typedef PolygonType::Pointer             PolygonPointerType;
  typedef PolygonType::ContinuousIndexType PolygonIndexType;
  typedef otb::ObjectList<PolygonType>     PolygonListType;
  typedef PolygonListType::Pointer         PolygonListPointerType;
  typedef unsigned long LabelPixelType;
  typedef otb::Image<LabelPixelType, 2>  LabeledImageType;

  typedef itk::ConnectedComponentImageFilter<CharImageType,LabeledImageType> ConnectedFilterType;
  typedef otb::PersistentVectorizationImageFilter<LabeledImageType,PolygonType> PersistentVectorizationFilterType;

  ConnectedFilterType::Pointer connectedFilter = ConnectedFilterType::New();
  connectedFilter->SetInput(l2bi->GetOutput());

  //Perform vectorization in a persistent way
    PersistentVectorizationFilterType::Pointer persistentVectorization = PersistentVectorizationFilterType::New();
    persistentVectorization->Reset();
    persistentVectorization->SetInput(connectedFilter->GetOutput());
    try
      {
      persistentVectorization->Update();
      }
    catch (itk::ExceptionObject& err)
      {
      std::cout << "\nExceptionObject caught !" << std::endl;
      std::cout << err << std::endl;
      }

	PolygonListPointerType OutputPolyList = persistentVectorization->GetPathList();
	//Display results
    std::cout << "nb objects found = " << OutputPolyList->Size() << std::endl;

	VectorDataType::Pointer outVectorData = VectorDataType::New(); 
    typedef VectorDataType::DataNodeType DataNodeType;
	typedef VectorDataType::DataTreeType            DataTreeType; 
    typedef itk::PreOrderTreeIterator<DataTreeType> TreeIteratorType;

	DataNodeType::Pointer document = DataNodeType::New(); 
    document->SetNodeType(otb::DOCUMENT); 
    document->SetNodeId("polygon"); 
    DataNodeType::Pointer folder = DataNodeType::New(); 
    folder->SetNodeType(otb::FOLDER); 
    DataNodeType::Pointer multiPolygon = DataNodeType::New(); 
    multiPolygon->SetNodeType(otb::FEATURE_MULTIPOLYGON);

	DataTreeType::Pointer tree = outVectorData->GetDataTree(); 
    DataNodeType::Pointer root = tree->GetRoot()->Get(); 
 
    tree->Add(document, root); 
    tree->Add(folder, document); 
    tree->Add(multiPolygon, folder);

	for (PolygonListType::Iterator pit = OutputPolyList->Begin(); 
       pit != OutputPolyList->End(); ++pit) 
    { 
    DataNodeType::Pointer newPolygon = DataNodeType::New(); 
    newPolygon->SetPolygonExteriorRing(pit.Get()); 
    tree->Add(newPolygon, multiPolygon); 
    }

	//outVectorData
   
	SetParameterOutputVectorData("out",outVectorData);    
    //SetParameterOutputImage<CharImageType>("out",l2bi->GetOutput());
  

  }



std::vector<CharImageType::IndexType> bresenham(CharImageType::IndexType start,
		CharImageType::IndexType end) {

	int x1 = start[0];
	int y1 = start[1];
	int x2 = end[0];
	int y2 = end[1];

	CharImageType::IndexType point;
	std::vector<CharImageType::IndexType> pointList;

	int delta_x(x2 - x1);
	// if x1 == x2, then it does not matter what we set here
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = std::abs(delta_x) << 1;

	int delta_y(y2 - y1);
	// if y1 == y2, then it does not matter what we set here
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = std::abs(delta_y) << 1;

	point[0] = x1;
	point[1] = y1;
	pointList.push_back(point);
	if (delta_x >= delta_y) {
		// error may go below zero
		int error(delta_y - (delta_x >> 1));

		while (x1 != x2) {
			if ((error >= 0) && (error || (ix > 0))) {
				error -= delta_x;
				y1 += iy;
			}
			// else do nothing

			error += delta_y;
			x1 += ix;

			point[0] = x1;
			point[1] = y1;
			pointList.push_back(point);

		}
	} else {
		// error may go below zero
		int error(delta_x - (delta_y >> 1));

		while (y1 != y2) {
			if ((error >= 0) && (error || (iy > 0))) {
				error -= delta_y;
				x1 += ix;
			}
			// else do nothing

			error += delta_x;
			y1 += iy;

			point[0] = x1;
			point[1] = y1;
			pointList.push_back(point);
		}
	}
	return pointList;
}

int verifier (int i, int j, int x, int y, itk::RescaleIntensityImageFilter<CharImageType,CharImageType>::Pointer image)

{
   unsigned int compte=0;

   //std::cout << " la valeur de i est  " << i << " la valeur de j est  " << j << std::endl;

  for (int i1 = -1; i1 < 2; i1++)
	   {
        for (int j1 = -1; j1 < 2; j1++)
	      {  
			
		      if ( i+i1 > -1 && i+i1 < x && j+j1 > -1 && j+j1 < y )
			  {  
                 if (i1 != 0 || j1 !=0)
				 {
			     CharImageType::IndexType pixelIndex; 

	             pixelIndex[0] = i+i1;   // x position 
                 pixelIndex[1] = j+j1;   // y position
	             CharImageType::PixelType pixelValue = image->GetOutput()->GetPixel(pixelIndex);

				 if (pixelValue == 255)
                 compte++; 
				 }
			  }
		   }  
	   }


	   if(compte==1)
       return 1;
	   return 0;
}


std::vector<std::vector<CharImageType::IndexType>> Recherche (itk::RescaleIntensityImageFilter<CharImageType,CharImageType>::Pointer image)
{   
	std::vector<std::vector<CharImageType::IndexType>> V;

	unsigned int compte=0;
	
	unsigned int x = image->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
	unsigned int y = image->GetOutput()->GetLargestPossibleRegion().GetSize()[1];


	for (unsigned int i = 0 ; i < x ; i++ )
	{	
		for (unsigned int j = 0 ; j < y ; j++ )

		{
	        CharImageType::IndexType pixelIndex; 

	        pixelIndex[0] = i;   // x position 
            pixelIndex[1] = j;   // y position
	        CharImageType::PixelType pixelValue = image->GetOutput()->GetPixel(pixelIndex);

	   if (pixelValue == 255 )
	     { 

	       if ( verifier(i,j,x,y,image) == 1  )	  
	       {	 
              for (int k=-100; k < 100 ; k++)
			  {
			     for (int l =-50; l < 50; l++ )
				 {  
				     if ( k+i>=0 && k+i<x && l+j>=0 && l+j<y )
					 
					{ 

                    CharImageType::IndexType pixelIndex; 
                    pixelIndex[0] = i+k;   // x position 
                    pixelIndex[1] = j+l;   // y position
	                CharImageType::PixelType pixelValue = image->GetOutput()->GetPixel(pixelIndex);

				    if ( (pixelValue == 255) && (verifier(i+k,j+l,x,y,image) == 1)  )	  
	                {
						

						CharImageType::IndexType pixelIndex1; 
                        pixelIndex1[0] = i;   // x position 
                        pixelIndex1[1] = j;   // y position

						CharImageType::IndexType pixelIndex2; 
                        pixelIndex2[0] = i+k;   // x position 
                        pixelIndex2[1] = j+l;   // y position

						std::vector<CharImageType::IndexType> pointList;
						pointList = bresenham(pixelIndex1, pixelIndex2);
						V.push_back(pointList);
				
					}

					}
				 }
			  }

			compte++;
	       }
	     }
		
		}}

	return V;
}


     };
}
   
}
OTB_APPLICATION_EXPORT(otb::Wrapper::SarRiverDetection);