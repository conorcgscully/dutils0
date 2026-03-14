import boto3
import pytest
import responses
from moto import mock_aws
from responses.registries import OrderedRegistry


@pytest.fixture
def aws():
    with mock_aws():
        yield


@pytest.fixture
def s3(aws):
    yield boto3.resource("s3", region_name="eu-central-1")


@pytest.fixture
def s3_bucket(s3):
    yield s3.create_bucket(
        Bucket="example-bucket",
        CreateBucketConfiguration={"LocationConstraint": "eu-central-1"},
    )


@pytest.fixture
def mocked_web():
    with responses.RequestsMock(registry=OrderedRegistry) as rsps:
        yield rsps


@pytest.fixture
def localfs(fs):
    # Wrap pyfakefs's `fs` fixture to avoid confusion with bindgen.fs
    return fs
