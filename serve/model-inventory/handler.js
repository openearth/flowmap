'use strict';

var AWS = require('aws-sdk');
var s3 = new AWS.S3();
var batch = new AWS.Batch();

module.exports.submit = (event, context, callback) => {

  const params = {
    jobDefinition: "matroos_flowmap",
    jobName: "matroos_flowmap",
    jobQueue: "flowmap-queue"
  };
  batch.submitJob(params, function(err, data) {
    if (err) {
      // an error occurred
      console.log(err, err.stack);
    } else {
      // successful response
      console.log(data);
    }
    const response = {
      statusCode: 200,
      body: JSON.stringify([
        data
      ])
    };
    callback(null, response);
  });


};
