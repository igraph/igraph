/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_UTIL_H
#define IGRAPHPP_UTIL_H

namespace igraph {

namespace noncopyable_ {

class noncopyable {
protected:
    noncopyable() {}
    ~noncopyable() {}

private:
    noncopyable(const noncopyable &);
    const noncopyable& operator=(const noncopyable&);
};

}

typedef noncopyable_::noncopyable noncopyable;

}         // end of namespaces

#endif    // IGRAPHPP_UTIL_H

