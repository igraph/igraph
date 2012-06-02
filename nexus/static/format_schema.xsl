<?xml version="1.0" encoding="UTF-8" ?>
<xsd:schema xmlns:xsd="http://www.w3.org/2001/XMLSchema">

  <xsd:annotation>
    <xsd:documentation xml:lang="en">
      Nexus network repository, data format(s)
      http://nexus.igraph.org 
    </xsd:documentation>
  </xsd:annotation>

  <xsd:element name="dataformats" type="dataformatsType" />
  
  <xsd:complexType name="dataformatsType">
    <xsd:sequence>
      <xsd:element name="dataformat" type="dataformatType" minOccurs="0"
		   maxOccurs="unbounded" />
    </xsd:sequence>
  </xsd:complexType>
  
  <xsd:complexType name="dataformatType">
    <xsd:sequence>
      <xsd:element name="name" type="xsd:string" />
      <xsd:element name="shortdescription" type="xsd:string" />
      <xsd:element name="description" type="xsd:string" />
      <xsd:element name="url" type="xsd:anyURI" />
    </xsd:sequence>
  </xsd:complexType>

</xsd:schema>
