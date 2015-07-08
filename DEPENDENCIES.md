#How to install Python packages

Because YAP uses its own copy of Python, you need to make sure that whatever additional
packages are installed, they are installed into that copy. Right now YAP does not have any
automated installation procedure, so the dependencies must be installed manually.

The proper command that uses easy_install would look like this (installing package
`argh` as an example and assuming `YAP_DEPS` is
already defined environment variable from YAP.rc and that the current Python version is
`2.7.4`:

```
$YAP_DEPS/python $YAP_DEPS/python-2.7.4/bin/easy_install -U argh
```

#Python package dependencies that must be installed

- argh
- numpy
- biopython
- setuptools

If suspecting that the list above in incomplete, check the contents of
`$YAP_DEPS/python-2.7.4/lib/python2.7/site-packages/

After the packages were installed, the Python directory is relocatable (meaning
that you can copy it into a different place and its Python executable will still
find the packages).

