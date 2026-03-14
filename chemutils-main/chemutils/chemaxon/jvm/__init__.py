"""
Provides access to a JVM running Chemaxon.

This uses the `py4j` library, which can start up a Java Virtual Machine which has access
to the Chemaxon Java code. The JVM is accessed by a call to `get_jvm()`, but alternatively
you can use the import magic of this package instead.

The magic means that instead of the normal use of `py4j`, where you would have to do:

```
from chemutils.chemaxon.jvm
logDPlugin = get_jvm().chemaxon.marvin.calculations.logDPlugin
```

you can instead do

```
from chemutils.chemaxon.jvm.chemaxon.marvin.calculations import logDPlugin
```
"""

import operator
import os
import sys
from collections.abc import Sequence
from functools import cache
from importlib.abc import Loader, MetaPathFinder
from importlib.machinery import ModuleSpec
from importlib.util import spec_from_loader
from types import ModuleType
from typing import Any, TypeGuard, TypeVar

import fsutils as fs
from py4j.java_gateway import JavaGateway, JavaObject, JVMView, is_instance_of

_T = TypeVar("_T")


class ChemaxonNotInstalledError(RuntimeError):
    def __init__(self) -> None:
        super().__init__(self, "Chemaxon is not installed.")


HAS_CHEMAXON = os.getenv("CHEMAXON_HOME") is not None


@cache
def get_java_gateway() -> JavaGateway:
    """Get a JVM with Chemaxon provided, caching so only one exists for the lifetime of the Python program."""
    if not HAS_CHEMAXON:
        raise ChemaxonNotInstalledError()

    existing_license_path = f"{os.getenv('CHEMAXON_HOME')}/license.cxl"
    if fs.file_exists(existing_license_path):
        license_path = existing_license_path
    else:
        license_path = fs.path.parse_path_as_absolute(
            "~/.cache/charm/chemaxon/license.cxl"
        ).geturl()
        try:
            fs.copy_file(
                src="s3://charmtx-datalake/licenses/chemaxon/license.cxl", dest=license_path
            )
        except OSError as e:
            raise RuntimeError(
                "Failed to copy Chemaxon license from S3, and no existing license was found."
            ) from e

    # Should not be None, because of the check in `chemutils.chemaxon` at import time
    chemaxon_path = os.getenv("CHEMAXON_HOME")
    assert chemaxon_path is not None

    return JavaGateway.launch_gateway(
        # Additional location to look for `.jar` files
        classpath=f"{chemaxon_path}/lib/*",
        # Provide the license
        javaopts=[f"-Dchemaxon.license.url={license_path}"],
        # Close the JVM when the Python thread is closed
        die_on_exit=True,
    )


@cache
def get_jvm() -> JVMView:
    """Get a JVM with Chemaxon provided, caching so only one exists for the lifetime of the Python program."""
    if not HAS_CHEMAXON:
        raise ChemaxonNotInstalledError()
    return get_java_gateway().jvm


def isjavainstance(obj: Any, clss: type[_T]) -> TypeGuard[_T]:
    """Is a given object an instance of a Java class."""
    if not isinstance(obj, JavaObject):
        return False
    if not HAS_CHEMAXON:
        raise ChemaxonNotInstalledError()
    return is_instance_of(get_java_gateway(), obj, clss)  # type: ignore


class ChemaxonMissingModule:
    def __init__(self) -> None:
        pass

    def __getattr__(self, name: str) -> Any:
        if name.startswith("__"):
            return super().__getattr__(self, name)  # type: ignore
        return ChemaxonMissingObject()


class ChemaxonMissingObject:
    def __getattr__(self, name: str) -> Any:
        raise ChemaxonNotInstalledError()

    def __call__(self) -> Any:
        raise ChemaxonNotInstalledError()


class JVMLoader(Loader):
    def __init__(self, attr: str):
        self.attr = attr

    def create_module(self, spec: ModuleSpec) -> Any:
        # Access the attribute as `jvm.a.b.c`
        # `operator.attrgetter` is like `getattr` but supports nested attributes
        obj = operator.attrgetter(self.attr)(get_jvm()) if HAS_CHEMAXON else ChemaxonMissingModule()

        # Any object with a `__path__` is treated by Python as a package
        obj.__path__ = __name__ + "." + self.attr  # type: ignore

        return obj

    def exec_module(self, module: ModuleType) -> None:
        # This is required to be an empty function
        pass


# Finder which can override how Python finds specific modules
# This overrides importing anything from a subpackage of `chemutils.chemaxon.jvm` as a
# request to import from a JVM itself
class JVMMetaPathFinder(MetaPathFinder):
    def find_spec(
        self, fullname: str, path: Sequence[str] | None = None, target: ModuleType | None = None
    ) -> Any:
        # Only handle things relative to this package
        if fullname.startswith(f"{__name__}."):
            return spec_from_loader(fullname, JVMLoader(fullname.removeprefix(f"{__name__}.")))


# Append the JVMMetaPathFinder to the system's `meta_path`, so that it is called when finding other packages.
sys.meta_path = [JVMMetaPathFinder(), *sys.meta_path]
