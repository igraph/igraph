REPLACE ----- remove the " * " prefix first -----------------*- mode:python -*-
^[ ]\*[ ]
WITH --------------------------------------------------------------------------
REPLACE ----- remove the " *" lines -------------------------------------------
^[ ]\*\s*\n
WITH --------------------------------------------------------------------------
\n
REPLACE ----- for the template functions --------------------------------------

FUNCTION\(
(?P<base>[^, \)]*)\s*,\s*
(?P<suffix>[^\)]*)
\)\s*

WITH

\g<base>_\g<suffix>

REPLACE ----- template type ---------------------------------------------------

TYPE\(
(?P<type>[^\)]*)
\)

WITH

\g<type>_t

REPLACE ----- template base type, we cowardly assume real number --------------

BASE

WITH 

igraph_real_t

REPLACE ----- function object, extract its signature --------------------------

(?P<before>\A.*?)                # head of the comment
\\function\s+                    # \function keyword
(?P<name>(?P<pre>(igraph_)|(IGRAPH_)|())(?P<tail>\w+)) # the keyword, remove igraph_ prefix
[\s]*(?P<brief>[^\n]*?)\n        # brief description
(?P<after>.*?)\*\/               # tail of the comment
\s*
(IGRAPH_EXPORT )?                # strip IGRAPH_EXPORT from prototype
(?P<def>.*?\))                   # function head
(?=(\s*;)|(\s*\{))               # prototype ends with ; function head with {
.*\Z                             # and the remainder

WITH --------------------------------------------------------------------------

<section id="\g<name>">
<title><function>\g<name></function> &mdash; \g<brief></title>
<indexterm><primary>\g<tail></primary></indexterm>
<para>
<informalexample><programlisting>
\g<def>;
</programlisting></informalexample>
</para>
<para>
\g<before>
\g<after>
</para>
</section>

REPLACE ----- <paramdef> for functions (not used currently) -------------------

<paramdef>(?P<params>[^<]*)</paramdef>\n

RUN ---------------------------------------------------------------------------

dr_params=string.split(matched.group("params"), ',')
dr_out=""
for dr_i in dr_params:
    dr_i=string.strip(dr_i)
    if dr_i=="...":
        dr_out=dr_out+"<varargs/>"
    else:
        dr_words=re.match(r"([\w\*\&\s]+)(\b\w+)$", dr_i).groups()
        dr_out=dr_out+"<paramdef>"+dr_words[0]+"<parameter>"+dr_words[1]+ \
                "</parameter></paramdef>\n"
actch=actch[0:matched.start()]+dr_out+actch[matched.end():]

REPLACE ----- function parameter descriptions, head ---------------------------

(?P<before>\A.*?)               # head of the comment
\\param\b                       # first \param commant

WITH --------------------------------------------------------------------------

\g<before></para>
<formalpara><title>Arguments:</title><para>
<variablelist role="params">
\\param

REPLACE ----- function parameter descriptions, tail ---------------------------

# the end of the params is either an empty line after the last \param
# command or a \return or \sa statement (others might be added later)
# or the end of the comment

\\param\b                         # the last \param command
(?P<paramtext>.*?)                # the text of the \param command
(?P<endmark>                      # this marks the end of the \param text
 (\\return\b)|(\\sa\b)|           # it is either a \return or \sa or
 (\n\s*?\n)|                      # (at least) one empty line or
 (\*\/))                          # the end of the comment
(?P<after>.*?\Z)                  # remaining part

WITH

\\param\g<paramtext></variablelist></para></formalpara><para>
\g<endmark>\g<after>

REPLACE ----- function parameter descriptions ---------------------------------

\\param\b\s*                      # \param command
(?P<paramname>(\w+)|(...))\s+     # name of the parameter
(?P<paramtext>.*?)                # text of the \param command
(?=(\\param)|(</variablelist>)|
 (\n\s*\n))


WITH --------------------------------------------------------------------------

  <varlistentry><term><parameter>\g<paramname></parameter>:</term>
  <listitem><para>
  \g<paramtext></para></listitem></varlistentry>  

REPLACE ----- \return command -------------------------------------------------

# a return statement ends with an empty line or the end of the comment
\\return\b\s*                     # \return command
(?P<text>.*?)                     # the text
(?=(\n\s*?\n)|                    # empty line or 
 (\*\/)|                          # the end of the comment or
 (\\sa\b))                        # \sa command

WITH ----------------------------------------------------------------------TODO

</para><formalpara><title>Returns:</title><para><variablelist>
  <varlistentry><term><parameter></parameter></term>
  <listitem><para>
  \g<text>
  </para></listitem></varlistentry>
</variablelist></para></formalpara><para>

REPLACE ----- variables -------------------------------------------------------

(?P<before>\A.*?)                 # head of the comment
\\var\s+                          # \var keyword + argument
(?P<name>(?P<pre>(igraph_)|(IGRAPH_)|())(?P<tail>\w+))
[\s]*(?P<brief>[^\n]*?)\n         # brief description
(?P<after>.*?)\*\/                # tail of the comment
\s*(?P<def>[^;]*;)                # the definition of the variable
.*\Z                              # and the remainder

WITH --------------------------------------------------------------------------

<section id="\g<name>"><title><function>\g<name></function> &mdash; \g<brief></title>
<indexterm><primary>\g<tail></primary></indexterm>
<para>
<programlisting>
\g<def>
</programlisting>
</para><para>
\g<before>\g<after>
</para>
</section>

REPLACE ----- \define ---------------------------------------------------------

(?P<before>\A.*?)                 # head of the comment
\\define\s+                       # \define command
(?P<name>(?P<pre>(igraph_)|(IGRAPH_)|())(?P<tail>\w+))
[\s]*(?P<brief>[^\n]*?)\n         # brief description
(?P<after>.*?)\*\/                # tail of the comment
\s*                               # whitespace
(?P<def>\#define\s+[\w0-9,]+\s*   # macro name
(\([\w0-9, ]+\))?)                # macro args (optional)
.*\Z                              # drop the remainder

WITH --------------------------------------------------------------------------

<section id="\g<name>"><title><function>\g<name></function> &mdash; \g<brief></title>
<indexterm><primary>\g<tail></primary></indexterm>
<para>
<programlisting>
\g<def>
</programlisting>
</para><para>
\g<before>\g<after>
</para>
</section>

REPLACE ----- \section without title ------------------------------------------

(?P<before>\A.*?)                 # head of the comment
\\section\s+(?P<name>\w+)\s*$     # \section + argument
(?P<after>.*?)\*\/                # tail of the comment
.*\Z                              # and the remainder, this is dropped

WITH

\g<before>
\g<after>

REPLACE ----- \section with title ---------------------------------------------

(?P<before>\A.*?)                 # head of the comment
\\section\s+(?P<name>\w+)         # \section + argument
(?P<title>.*?)                    # section title
\n\s*?\n                          # empty line
(?P<after>.*?)\*\/                # tail of the comment
.*\Z                              # and the remainder, this is dropped

WITH

<title>\g<title></title>
\g<before>
\g<after>

REPLACE ----- \section with title ---------------------------------------------

(?P<before>\A.*?)                 # head of the comment
\\section\s+(?P<name>\w+)         # \section + argument
(?P<title>.*?)\s*\*\/             # section title
.*\Z                              # and the remainder, this is dropped

WITH

<title>\g<title></title>
\g<before>

REPLACE ----- an enumeration typedef ------------------------------------------

(?P<before>\A.*?)                 # head of the comment
\\typedef\s+                      # \typedef command
(?P<name>(?P<pre>(igraph_)|(IGRAPH_)|())(?P<tail>\w+))
[\s]*(?P<brief>[^\n]*?)\n         # brief description
(?P<after>.*?)                    # tail of the comment
 \*\/\s*                          # closing the comment
(?P<def>typedef\s*enum\s*\{       # typedef enum
 [^\}]*\}\s*\w+\s*;)                  # rest of the definition
.*\Z

WITH --------------------------------------------------------------------------

<section id="\g<name>"><title><function>\g<name></function> &mdash; \g<brief></title>
<indexterm><primary>\g<tail></primary></indexterm>
<para>
<programlisting>
\g<def>
</programlisting>
</para>
<para>
\g<before>\g<after>
</para>
</section>

REPLACE ----- enumeration value descriptions, head ----------------------------

(?P<before>\A.*?)               # head of the comment
\\enumval\b                     # first \param commant

WITH --------------------------------------------------------------------------

\g<before></para>
<formalpara><title>Values:</title><para>
<variablelist role="params">
\\enumval

REPLACE ----- enumeration value descriptions, tail ----------------------------

\\enumval\b                       # the last \enumval command
(?P<paramtext>.*?)                # the text of the \enumval command
(?P<endmark>                      # this marks the end of the \enumval text
 (\\return\b)|(\\sa\b)|           # it is either a \return or \sa or
 (\n\s*?\n)|                      # (at least) one empty line or
 (\*\/))                          # the end of the comment
(?P<after>.*?\Z)                  # remaining part

WITH

\\enumval\g<paramtext></variablelist></para></formalpara><para>
\g<endmark>\g<after>

REPLACE ----- enumeration value descriptions ----------------------------------

\\enumval\b\s*                    # \enumval command
(?P<paramname>(\w+)|(...))\s+     # name of the parameter
(?P<paramtext>.*?)                # text of the \enumval command
(?=(\\enumval)|(</variablelist>)|
 (\n\s*\n))

WITH --------------------------------------------------------------------------

  <varlistentry><term><constant>\g<paramname></constant>:</term>
  <listitem><para>
  \g<paramtext></para></listitem></varlistentry>  

REPLACE ----- \struct ---------------------------------------------------------

(?P<before>\A.*?)                 # head of the comment
\\struct\s+                       # \struct command
(?P<name>(?P<pre>(igraph_)|(IGRAPH_)|())(?P<tail>[\w_]+))
[\s]*(?P<brief>[^\n]*?)(?=\n)     # brief description
(?P<after>.*?)                    # tail of the command
\*\/\s*                           # closing the comment
(?P<def>typedef \s*struct\s*\w+\s*\{
 .*\}\s*\w+\s*;)
.*\Z

WITH --------------------------------------------------------------------------

<section id="\g<name>"><title><function>\g<name></function> &mdash; \g<brief></title>
<indexterm><primary>\g<tail></primary></indexterm>
<para>
<programlisting>
\g<def>
</programlisting>
</para>
<para>
\g<before>\g<after>
</para>
</section>

REPLACE ----- structure member descriptions, one block ------------------------

^[\s]*\n
(?P<before2>.*?)                # empty line+text
(?P<members>\\member\b.*?)      # member commands
(?=                             # this marks the end of the \member text
 (\\return\b)|(\\sa\b)|         # it is either a \return or \sa or
 (^[\s]*\n)|                    # (at least) one empty line or
 (\*\/))                        # the end of the comment

WITH --------------------------------------------------------------------------

</para>
<para>\g<before2></para>
<formalpara><title>Values:</title>
<para><variablelist role="params">
\g<members>
</variablelist></para></formalpara><para>

REPLACE ----- structure member descriptions -----------------------------------

\\member\b\s*                    # \enumval command
(?P<paramname>(\w+)|(...))\s+     # name of the parameter
(?P<paramtext>.*?)                # text of the \enumval command
(?=(\\member)|(</variablelist>)|
 (\n\s*\n))

WITH --------------------------------------------------------------------------

  <varlistentry><term><constant>\g<paramname></constant>:</term>
  <listitem><para>
  \g<paramtext></para></listitem></varlistentry>

REPLACE ----- \typedef function -----------------------------------------------

(?P<before>\A.*?)                   # comment head
\\typedef\s+                      # \typedef command
(?P<name>(?P<pre>(igraph_)|(IGRAPH_)|())(?P<tail>\w+))
[\s]*(?P<brief>[^\n]*?)\n         # brief description
(?P<after>.*?)                    # comment tail
\*\/                              # end of comment block
\s*
(?P<src>typedef\s+[^;]*;)        # the typedef definition
.*\Z

WITH --------------------------------------------------------------------------

<section id="\g<name>"><title><function>\g<name></function> &mdash; \g<brief></title>
<indexterm><primary>\g<tail></primary></indexterm>
<para><programlisting>
\g<src>
</programlisting></para>
<para>
\g<before>\g<after>
</para>
</section>

REPLACE ----- ignore doxygen \ingroup command ---------------------------------

\\ingroup\s+\w+

WITH --------------------------------------------------------------------------

REPLACE ----- ignore doxygen \defgroup command --------------------------------

\\defgroup\s+\w+

WITH --------------------------------------------------------------------------

REPLACE ----- add the contents of \brief to the description -------------------

\\brief\b

WITH --------------------------------------------------------------------------

REPLACE ----- \varname command ------------------------------------------------

\\varname\b\s*
(?P<var>\w+\b)

WITH

<varname>\g<var></varname>

REPLACE ----- references, \ref command ----------------------------------------

\\ref\b\s*
(?P<what>\w+)(?P<paren>([\(][\)])?)

WITH --------------------------------------------------------------------------

<link linkend="\g<what>"><function>\g<what>\g<paren></function></link>

REPLACE ----- \sa command -----------------------------------------------------

\\sa\b
\s*
(?P<text>.*?)
(?=(\n\s*?\n)|(\*\/))

WITH ----------------------------------------------------------------------TODO

</para><formalpara><title>See also:</title><para><variablelist>
  <varlistentry><term><parameter></parameter></term>
  <listitem><para>
  \g<text>
  </para></listitem></varlistentry>
</variablelist></para></formalpara><para>

REPLACE ----- \em command -----------------------------------------------------

\\em\b
\s*
(?P<text>[^\s]+)

WITH

<emphasis>\g<text></emphasis>

REPLACE ----- \emb command ----------------------------------------------------

\\emb\b

WITH

<emphasis>

REPLACE ----- \eme command ----------------------------------------------------

\\eme\b

WITH

</emphasis>

REPLACE ----- \verbatim -------------------------------------------------------

\\verbatim\b

WITH

<informalexample><programlisting>

REPLACE ----- \endverbatim ----------------------------------------------------

\\endverbatim\b

WITH

</programlisting></informalexample>

REPLACE ----- \clist ----------------------------------------------------------

\\clist\b

WITH

<variablelist>

REPLACE ----- \cli ------------------------------------------------------------

\\cli\s+(?P<term>.*?)$
(?P<text>.*?)
(?=(\\cli)|(\\endclist))

WITH --------------------------------------------------------------------------

<varlistentry><term><constant>\g<term></constant></term>
<listitem><para>
\g<text>
</para></listitem></varlistentry>

REPLACE ----- \endclist -------------------------------------------------------

\\endclist\b

WITH

</variablelist>

REPLACE ----- \olist ----------------------------------------------------------

\\olist\b

WITH

<orderedlist>

REPLACE ----- \oli ------------------------------------------------------------

\\oli\s+(?P<text>.*?)
(?=(\\oli)|(\\endolist))

WITH

<listitem><para>
\g<text>
</para></listitem>

REPLACE ----- \endolist -------------------------------------------------------

\\endolist\b

WITH

</orderedlist>

REPLACE ----- \ilist ----------------------------------------------------------

\\ilist\b

WITH

<itemizedlist>

REPLACE ----- \ili ------------------------------------------------------------

\\ili\s+(?P<text>.*?)
(?=(\\ili)|(\\endilist))

WITH

<listitem><para>
\g<text>
</para></listitem>

REPLACE ----- \endilist -------------------------------------------------------

\\endilist\b

WITH

</itemizedlist>

REPLACE ----- doxygen \c command is for <constant> ----------------------------

\\c\s+(?P<word>[\w\-^\']+)\b

WITH

<constant>\g<word></constant>

REPLACE ----- doxygen \p command is for <parameter> ---------------------------

\\p\s+(?P<word>\w+)\b

WITH

<parameter>\g<word></parameter>

REPLACE ----- doxygen \type command is for <type> -----------------------------

\\type\s+(?P<word>\w+)\b

WITH

<type>\g<word></type>

REPLACE ----- doxygen \a command is for <command> -----------------------------

\\a\s+(?P<word>\w+)\b

WITH

<command>\g<word></command>

REPLACE ----- doxygen \quote command is for <quote> ---------------------------

\\quote\s+

WITH

<quote>

REPLACE ----- doxygen \endquote command is for </quote> -----------------------

\s*\\endquote\b

WITH

</quote>

REPLACE ----- replace <code> with <literal> -----------------------------------

<(?P<c>/?)code>

WITH --------------------------------------------------------------------------

<\g<c>literal> 

REPLACE ----- add http:// and https:// links ----------------------------------

(?P<link>https?:\/\/[-\+=&;%@./~()'\w_]*[-\+=&;%@/~'\w_])

WITH --------------------------------------------------------------------------

<ulink url="\g<link>">\g<link></ulink>

REPLACE ----- blockquote ------------------------------------------------------

\\blockquote

WITH --------------------------------------------------------------------------

<blockquote>

REPLACE ----- blockquote ------------------------------------------------------

\\endblockquote

WITH --------------------------------------------------------------------------

</blockquote>

REPLACE ----- example file  ---------------------------------------------------

\\example\b\s*
(?P<filename>[^\n]*?)\n

WITH --------------------------------------------------------------------------

<example role="sourcefile">
  <title> File <code>\g<filename></code></title>
  <xi:include  href="../\g<filename>.xml"
      xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <para></para>
</example>

REPLACE ----- \deprecated-by --------------------------------------------------

\\deprecated-by\b\s*
(?P<replacement>[^ \n]+)\s*
(?P<version>[^\n]+)\n

WITH --------------------------------------------------------------------------

</para>
<warning>
<para>Deprecated since version \g<version>. Please do not use this function in new
code; use <link linkend="\g<replacement>"><function>\g<replacement>()</function></link>
instead.</para>
</warning>
<para>

REPLACE ----- \deprecated -----------------------------------------------------

\\deprecated\b\s*
(?P<version>[^\n]*?)\n

WITH --------------------------------------------------------------------------

</para>
<warning>
<para>Deprecated since version \g<version>. Please do not use this function in new
code.</para>
</warning>
<para>

REPLACE ----- \experimental ---------------------------------------------------

\\experimental\b\s*\n

WITH --------------------------------------------------------------------------

</para>
<warning>
<para>This function is experimental and its signature is not considered final yet.
We reserve the right to change the function signature without changing the
major version of igraph. Use it at your own risk.</para>
</warning>
<para>

