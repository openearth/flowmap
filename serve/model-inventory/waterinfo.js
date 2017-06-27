'use strict';
var fetch = require('node-fetch');
var proj4 = require('proj4');
var reproject = require('reproject').reproject;
var _ = require('lodash');

// define rd
proj4.defs['EPSG:28992'] = new proj4.Proj('+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs');
proj4.defs['EPSG:25831'] = new proj4.Proj('+proj=utm +zone=31 +ellps=GRS80 +units=m +no_defs');

const headers = {
    'Access-Control-Allow-Origin' : '*', // Required for CORS support to work
    'Access-Control-Allow-Credentials' : true // Required for cookies, authorization headers with HTTPS
};

function points(event, context, callback) {
    // urls
    const pointsUrl = 'http://waterinfo.rws.nl/api/point/latestmeasurements?parameterIds=Waterhoogte___20Oppervlaktewater___20t.o.v.___20Normaal___20Amsterdams___20Peil___20in___20cm';

    // projections
    let rd = proj4('EPSG:28992');
    let utm = proj4('EPSG:25831');
    let wgs84 = proj4('EPSG:4326');

    // get data
    fetch(pointsUrl)
        .then(function(resp) {
            return resp.json();
        })
        .then(function(json) {
            let projected = reproject(json, utm, wgs84);

            const featured = [
                // TODO: fix these
                'IJmuiden Buitenhaven', 'Delfzijl', 'Den Helder',
                'Vlissingen'
            ];
            _.each(projected.features, function(feature) {
                // update featured items
                let name = feature.properties.name;
                feature.properties.featured = _.includes(featured, name);
            });
            const response = {
                statusCode: 200,
                headers: headers,
                body: JSON.stringify(projected)
            };
            callback(null, response);

        })
        .catch(function(err) {
            console.warn(err);
        });

}

function details(event, context, callback) {
    let location = _.get(
        event,
        // get the pat parameter
        'pathParameters.location',
        // or location if path parameters is unavailable
        event.location
    );

    // fill in url
    let detailsUrl = `http://waterinfo.rws.nl/api/details/chart?mapType=waterhoogte-t-o-v-nap&values=-240,48&locationCode=${location}`;

    // fetch it and return it
    fetch(detailsUrl)
        .then(resp => {
            return resp.json();
        })

        .then(function(json) {
            // convert from cm to m
            _.each(json.limits, function(limit) {
                if (_.isNumber(limit.from))  {
                    limit.from = limit.from / 100.0;
                }
                if (_.isNumber(limit.to))  {
                    limit.to = limit.to / 100.0;
                }
                limit.label = limit.label.replace(/\bcm\b/gi, 'm');
            });
            // convert series from cm to m
            _.each(json.series, function(serie) {
                _.each(serie.data, function(row) {
                    if (_.isNumber(row.value))  {
                        row.value = row.value / 100.0;
                    }
                });
                serie.name = serie.name.replace(/\bcm\b/g, 'm');
                serie.unit = serie.unit.replace(/\bcm\b/g, 'm');
            });
            if (_.has(json, 'extremesY')) {
                json.extremesY.min = json.extremesY.min/100.0;
                json.extremesY.max = json.extremesY.max/100.0;
            }

            const response = {
                statusCode: 200,
                headers: headers,
                body: JSON.stringify(json)
            };
            callback(null, response);
        })
        .catch(function(err) {
            console.warn(err);
        });
}

module.exports = {
    points,
    details
};
