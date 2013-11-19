"use strict"

var express = require('express'),
  hbs = require('hbs');

var app = express();
app.set('view engine', 'html');
app.engine('html', require('hbs').__express);
app.set('views', __dirname + '/views');
app.use(express.static(__dirname + '/public'));
app.use(express.logger());

app.get('/', function(req, res){
  res.render('index_template', {});
});

app.get('/matrix', function(req, res){
    res.render('matrix_template', {})
});

var port = process.env.PORT || 5000;


app.listen(port, function() {
  console.log("listening on " + port);
})
