<?xml version="1.0" encoding="UTF-8" ?>
<xsd:schema xmlns:xsd="http://www.w3.org/2001/XMLSchema">

  <xsd:annotation>
    <xsd:documentation xml:lang="en">
      Nexus network repository, licence information
      http://nexus.igraph.org 
    </xsd:documentation>
  </xsd:annotation>

  <xsd:element name="licences" type="licencesType" />

  <xsd:complexType name="licencesType">
    <xsd:sequence>
      <xsd:element name="licence" type="licenceType" minOccurs="0"
		   maxOccurs="unbounded" />
    </xsd:sequence>
  </xsd:complexType>
  
  <xsd:complexType name="licenceType">
    <xsd:sequence>
      <xsd:element name="id" type="xsd:positiveInteger" />
      <xsd:element name="name" type="xsd:string" />
      <xsd:element name="shortdescription" type="xsd:string" />
      <xsd:element name="url" type="xsd:anyURI" />
      <xsd:element name="text" type="xsd:string" minOccurs="0" />
    </xsd:sequence>
  </xsd:complexType>
  
</xsd:schema>

      