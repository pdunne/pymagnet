#!/bin/sh
SRC_PATH="./src/pymagnet"
REP_PATH="./reports"
COV_PATH="./reports/coverage"
JUNIT_PATH="./reports/junit"
FLAKE_PATH="./reports/flake8"
TEST_PATH="./tests"

poetry run pytest --cov pymagnet $TEST_PATH/ \
--cov-branch \
--cov-report term \
--cov-report html:$COV_PATH/html \
--cov-report xml:$COV_PATH/coverage.xml

poetry run genbadge coverage -i $COV_PATH/coverage.xml -o $REP_PATH/coverage-badge.svg

poetry run pytest --junitxml=$JUNIT_PATH/junit.xml 
poetry run genbadge tests -i $JUNIT_PATH/junit.xml  -o $REP_PATH/tests-badge.svg 

poetry run flake8 $SRC_PATH --ignore=E203,E266,E501,W503 --exit-zero --format=html --htmldir $FLAKE_PATH --statistics --tee --output-file $FLAKE_PATH/flake8stats.txt
poetry run genbadge flake8 -i $FLAKE_PATH/flake8stats.txt -o $REP_PATH/flake8-badge.svg