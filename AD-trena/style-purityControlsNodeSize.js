vizmap = [

   {selector: "node", css: {
      "shape": "ellipse",
      "text-valign":"center",
      "text-halign":"center",
      "content": "data(label)",
      "background-color": "#FFFFFF",
      "border-color":"black","border-width":"1px",
      "width": "mapData(degree, 0.0, 100.0, 20.0, 200.0)",
      "height":"mapData(degree, 0.0, 100.0, 20.0, 200.0)",
      "font-size":"32px"}},


   {selector: "node[type='info']", css: {
       "shape": "roundrectangle",
       "font-size": "72px",
       "width": "360px",
       "height": "120px",
       "border-width": "3px",
       "background-color": "beige"
       }},

   {selector: "node[type='footprint']", css: {
       "shape": "roundrectangle",
       "width": "40px",
       "height": "30px",
       "font-size": "14px",
       "background-color": "beige"
       }},

   {selector: "node[type='TF']", css: {
      "shape": "ellipse",
      "background-color": "mapData(beta, 0, 1.0, white, red)",
      "width": "mapData(purity, 0.0, 30.0, 20.0, 200.0)",
      "height":"mapData(purity, 0.0, 30.0, 20.0, 200.0)"
       }},

   {selector: "node[type='targetGene']", css: {
      "shape": "octagon",
      "background-color": "lightBlue",
      "width": 100,
      "height": 100
       }},

   {selector: "node[type='TF'][beta<=0]", css: {
      "shape": "ellipse",
      "background-color": "mapData(beta, -1.0, 0, green, white)"
      //"width": "mapData(purity, 0.0, 30.0, 20.0, 200.0)",
      //"height":"mapData(purity, 0.0, 30.0, 20.0, 200.0)"
       }},

   {selector: "node[type='tf']", css: {
       "shape": "rectangle"
       }},


   {selector: "node[score=1]", css: {
       "background-color": "#AAFFAA"
       }},

   {selector: "node[score=2]", css: {
       "background-color": "#DDDDFF"
       }},

   {selector: "node[score=3]", css: {
       "background-color": "#FF2222"
       }},

   {selector: 'edge', css: {
      "line-color": "rgb(200, 200, 200)",
      "target-arrow-shape": "triangle",
      "target-arrow-color": "rgb(0, 0, 0)",
      width: "1px"
      }},

    {selector: 'edge', css: {
       'line-color': 'black',
       'target-arrow-shape': 'triangle',
       'target-arrow-color': 'black',
       'curve-style': 'bezier'
        }},

   {selector: 'edge:selected', css: {
      "line-color": "red",
      "target-arrow-shape": "triangle",
      "target-arrow-color": "rgb(0, 0, 0)",
       width: "5px",
      'curve-style': 'bezier'
      }},

   {selector:"node:selected", css: {
       "text-valign":"center",
       "text-halign":"center",
       "border-color": "black",
       "content": "data(id)",
       "border-width": "3px",
       "overlay-opacity": 0.2,
       "overlay-color": "gray"
        }}


   ];
