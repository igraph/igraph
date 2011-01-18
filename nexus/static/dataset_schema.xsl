<?xml version="1.0" encoding="UTF-8" ?>
<xsd:schema xmlns:xsd="http://www.w3.org/2001/XMLSchema">

  <xsd:annotation>
    <xsd:documentation xml:lang="en">
      Nexus network repository, dataset information
      http://nexus.igraph.org 
    </xsd:documentation>
  </xsd:annotation>

  <xsd:element name="dataset" type="datasetType"/>

  <xsd:complexType name="datasetType">
    <xsd:sequence>
      <xsd:element name="id" type="xsd:positiveInteger" />
      <xsd:element name="sid" type="xsd:string" />
      <xsd:element name="name" type="xsd:string" />
      <xsd:element name="shortdescription" type="xsd:string" />
      <xsd:element name="description" type="xsd:string" />
      <xsd:element name="vertices" type="xsd:nonNegativeInteger" />
      <xsd:element name="edges" type="xsd:nonNegativeInteger" />
      <xsd:element name="tags" type="tagsType" minOccurs="0" />
      <xsd:element name="attributes" type="attributesType" />
      <xsd:element name="date" type="xsd:date" />
      <xsd:element name="licence" type="xsd:string" />
      <xsd:element name="licenceurl" type="xsd:anyURI" />
      <xsd:element name="papers" type="papersType" />
      <xsd:element name="formats" type="formatsType" />
    </xsd:sequence>
  </xsd:complexType>

  <xsd:complexType name="tagsType">
    <xsd:sequence>
      <xsd:element name="tag" type="xsd:string" minOccurs="0" 
		   maxOccurs="unbounded" />
    </xsd:sequence>
  </xsd:complexType>

  <xsd:complexType name="attributesType">
    <xsd:sequence>
      <xsd:element name="attribute" type="attributeType" minOccurs="0"
		   maxOccurs="unbounded" />
    </xsd:sequence>
  </xsd:complexType>

  <xsd:complexType name="attributeType">
    <xsd:sequence>
      <xsd:element name="name" type="xsd:string" />
      <xsd:element name="type" type="xsd:string" />
      <xsd:element name="datatype" type="xsd:string" />
      <xsd:element name="description" type="xsd:string" />
    </xsd:sequence>
  </xsd:complexType>

  <xsd:complexType name="papersType">
    <xsd:sequence>
      <xsd:element name="paper" type="xsd:string" minOccurs="0"
		   maxOccurs="unbounded" />
    </xsd:sequence>
  </xsd:complexType>  

  <xsd:complexType name="formatsType">
    <xsd:sequence>
      <xsd:element name="format" type="xsd:string" minOccurs="0"
		   maxOccurs="unbounded" />
    </xsd:sequence>
  </xsd:complexType>  

</xsd:schema>
