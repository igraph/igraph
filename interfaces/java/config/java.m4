# Detecting Java environment                  -*- Autoconf -*-

AC_DEFUN([AC_PROG_JAVA],
[
have_java=no

if test -z "${JAVA_HOME}"; then
  JAVA_PATH=${PATH}
else
  JAVA_PATH=${JAVA_HOME}:${JAVA_HOME}/jre/bin:${JAVA_HOME}/bin:${JAVA_HOME}/../bin:${PATH}
fi
AC_PATH_PROGS(JAVA, java, ,${JAVA_PATH})
AC_PATH_PROGS(JAR, jar, ,${JAVA_PATH})
AC_PATH_PROGS(JAVAC, javac, ,${JAVA_PATH})
AC_PATH_PROGS(JAVAH, javah, ,${JAVA_PATH})

cat >GetProp.java <<EOF
import java.util.StringTokenizer;
public class GetProp {
  public static void main(String[[]] args) {
    if (args==null || args.length==0) return;
    if (args[[0]].equals("--libs")) {
      String prefix="-L";
      String libpath=System.getProperty("java.library.path");
      String pathsep=System.getProperty("path.separator");
      if (pathsep==null || pathsep.length()==0) pathsep=":";
      if (args.length>1) prefix=args[[1]];
      StringTokenizer tok=new StringTokenizer(libpath,pathsep);
      while (tok.hasMoreTokens()) {
        System.out.print(prefix);
	System.out.print(tok.nextToken());
	if (tok.hasMoreTokens()) System.out.print(" ");
      }
      System.out.println("");
      return;
    }
    System.out.println(System.getProperty(args[[0]]));
  }
}
EOF

AC_CACHE_CHECK([whether Java compiler works], [igraph_java_works], [
  igraph_java_works=no
  if test -n "${JAVAC}"; then
    rm -f GetProp.class
    if "${JAVAC}" GetProp.java 2>&AS_MESSAGE_LOG_FD; then
      if test -e GetProp.class; then igraph_java_works=yes; fi
    fi
  fi
])

if test ${igraph_java_works} = yes; then
  AC_CACHE_CHECK([Java home], [igraph_java_home], [
    if test -z "${JAVA_HOME}"; then
      JAVA_HOME=`${JAVA} GetProp java.home`
    fi
    igraph_java_home="${JAVA_HOME}"
  ])
  JAVA_HOME="${igraph_java_home}"
  have_java=yes
fi

if test -n "${JAVA_HOME}"; then
  host_os=`${JAVA} GetProp os.name`
  case "${host_os}" in
    Mac*)
      JAVA_LIBS="-framework JavaVM"
      JAVA_CPPFLAGS="-dynamiclib -I${JAVA_HOME}/include"
      JAVA_LD_LIBRARY_PATH=
    ;;
    *)
      JAVA_LIBS="`${JAVA} GetProp -libs` -ljvm"
      JAVA_LD_LIBRARY_PATH="`${JAVA} GetProp java.library.path`"
      JAVA_CPPFLAGS="-I${JAVA_HOME}/include"
    ;;
  esac
  
  AC_CACHE_CHECK([whether JNI works], [igraph_jni_works],
  [igraph_jni_works=no
  igraph_save_CPPFLAGS="${CPPFLAGS}"
  igraph_save_LIBS="${LIBS}"
  LIBS="${JAVA_LIBS}"
  CPPFLAGS="${JAVA_CPPFLAGS}"
  AC_LINK_IFELSE([
#include <jni.h>
int main() { JNI_CreateJavaVM(0,0,0); return 0; }],
  [igraph_jni_works=yes],[igraph_jni_works=no])
  ])
fi

rm -rf GetProp.java GetProp.class

AC_SUBST(JAVA)
AC_SUBST(JAVAC)
AC_SUBST(JAVAH)
AC_SUBST(JAR)
AC_SUBST(JAVA_LD_LIBRARY_PATH)
AC_SUBST(JAVA_LIBS)
AC_SUBST(JAVA_CPPFLAGS)
])
