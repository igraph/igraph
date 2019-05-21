function getElementByClass(element, className) {
    tc=element.childNodes;
    for (var i=0; i<tc.length; i++) {
	if (tc[i].className == className) { return tc[i]; }
    }
    return null;
}

function toggle(target, show)
{
    exdiv=getElementByClass(target, "example");
    excdiv=getElementByClass(exdiv, "example-contents");

    if (excdiv.style.display != 'block')
	excdiv.style.display = 'block';
    else
	excdiv.style.display = 'none';
}
