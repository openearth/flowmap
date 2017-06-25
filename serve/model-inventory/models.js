'use strict';

var AWS = require('aws-sdk');
var s3 = new AWS.S3();
var _ = require('lodash');

function listObjects (prefix) {
  return new Promise(function(resolve, reject) {
    var params = {
      Bucket: "model-inventory",
      Prefix: prefix
    };
    s3.listObjects(params, function(err, data) {
      if (err) {
        reject(err);
      } // an error occurred
      else {
        resolve(data.Contents);
      }
    });
  });
};

function getObject(key) {
  return new Promise(function(resolve, reject) {
    var params = {
      Bucket: "model-inventory",
      Key: key
    };
    s3.getObject(params, function(err, data) {
      if (err) {
        reject(err);
      } // an error occurred
      else {
        var body = data.Body.toString('utf-8');
        var model = JSON.parse(body);
        resolve(model);
      }
    });
  });
};

module.exports.index = (event, context, callback) => {
  Promise.all([
    listObjects("meta/"),
    listObjects("movies/")
  ]).then((result) => {
    var metas = result[0];
    var movies = result[1];
    var models = [];
    var modelPromises = _.map(
      // only the non empty files
      _.filter(
        metas,
        (meta) => meta.Size > 0
      ),
      (meta) => {
        // get each object
        return getObject(meta.Key);
      }
    );

    Promise.all(modelPromises)
      .then(
        (models) => {
          const response = {
            statusCode: 200,
            body: JSON.stringify({
              "models": models
            })
          };
          callback(null, response);
        }
      );
  });
};
