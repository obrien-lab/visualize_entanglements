// Create NGL Stage object
var stage = new NGL.Stage( "viewport" );

// Handle window resizing
window.addEventListener( "resize", function( event ){
    stage.handleResize();
}, false );


// Code for example: interactive/simple-viewer
stage.setParameters({
  backgroundColor: "white"
})

function addElement (container, el) {
  Object.assign(el.style, {
    position: "absolute",
    zIndex: 10
  })
  container.appendChild(el)
}

function createElement (name, properties, style) {
  var el = document.createElement(name)
  Object.assign(el, properties)
  Object.assign(el.style, style)
  return el
}

function createSelect (options, properties, style) {
  var select = createElement("select", properties, style)
  options.forEach(function (d) {
    select.add(createElement("option", {
      value: d[ 0 ], text: d[ 1 ]
    }))
  })
  return select
}

function createFileButton (label, properties, style) {
  var input = createElement("input", Object.assign({
    type: "file"
  }, properties), { display: "none" })
  addElement(document.getElementById("cntrl_panel"), input)
  var button = createElement("input", {
    value: label,
    type: "button",
    onclick: function () { input.click() }
  }, style)
  return button
}

// recenter the structure in the scene
function reset_view() {
  var sele_1 = ""
  var sele_2 = ""
  stage.getRepresentationsByName("ligand")
    .setVisibility(ligandCheckbox.checked)
  if (ligandCheckbox.checked) {
    sele_2 += "(" + sele_dict["ligand"] + ") or "
  }
  stage.getRepresentationsByName("protein_2")
    .setVisibility(nativeCheckbox.checked)
  stage.getRepresentationsByName("circle_2")
    .setVisibility(nativeCheckbox.checked)
  stage.getRepresentationsByName("thread_2")
    .setVisibility(nativeCheckbox.checked)
  if (nativeCheckbox.checked) {
    sele_2 += "(" + sele_dict["protein_2"] + ")"
  }
  stage.getRepresentationsByName("protein_1")
    .setVisibility(stateCheckbox.checked)
  stage.getRepresentationsByName("ribosome")
    .setVisibility(stateCheckbox.checked)
  stage.getRepresentationsByName("tRNAs")
    .setVisibility(stateCheckbox.checked)
  if (stateCheckbox.checked) {
    sele_1 += "(" + sele_dict["protein_1"] + ") or (" + sele_dict["ribosome"] + ") or (" + sele_dict["tRNAs"] + ") or "
  }
  stage.getRepresentationsByName("circle_1")
    .setVisibility(entanglementCheckbox.checked)
  stage.getRepresentationsByName("thread_1")
    .setVisibility(entanglementCheckbox.checked)
  stage.getRepresentationsByName("NC_atoms")
    .setVisibility(entanglementCheckbox.checked)
  stage.getRepresentationsByName("native contact")
    .setVisibility(entanglementCheckbox.checked)
  if (entanglementCheckbox.checked) {
    sele_1 += "(" + sele_dict["circle_1"] + ") or (" + sele_dict["thread_1"] + ")"
  }

  if (sele_1.length > 0 && sele_1[sele_1.length-1] == " ") {
    sele_1 = sele_1.substr(0, sele_1.length-4)
  }
  if (sele_2.length > 0 && sele_2[sele_2.length-1] == " ") {
    sele_2 = sele_2.substr(0, sele_2.length-4)
  }

  sele_1 = sele_1.trim()
  sele_2 = sele_2.trim()
  
  var box_1 = new THREE.Box3()
  var box_2 = new THREE.Box3()

  if (sele_1.length > 0) {
    box_1 = struct[0].structure.getBoundingBox(new NGL.Selection(sele_1))
  }

  if (sele_2.length > 0) {
    box_2 = struct[1].structure.getBoundingBox(new NGL.Selection(sele_2))
  }

  var b_min = box_1.min.clone()
  b_min = b_min.min(box_2.min)
  var b_max = box_1.max.clone()
  b_max = b_max.max(box_2.max)
  var new_box = new THREE.Box3(b_min, b_max)

  if (! new_box.isEmpty()) {
    stage.animationControls.zoomMove(
      new_box.getCenter(),
      stage.getZoomForBox(new_box),
      1000
    )
  }
}

// calculate gN and gC
function calc_g(g_list, M, resid_1, resid_2) {
  var terminal_cutoff = 5
  var range = [[terminal_cutoff, (resid_1-1)-4], [(resid_2-1)+4, M.length-terminal_cutoff]]
  var i;
  for (i = 0; i < 2; i++) {
    var j1;
    var g = 0;
    for (j1 = resid_1-1; j1 < resid_2-1; j1++) {
      var j2;
      for (j2 = range[i][0]; j2 < range[i][1]; j2++) {
        g += M[j1][j2];
      }
    }
    g = Math.round(g/4/3.14)
    g_list.push(g)
  }
}

// Generate the topology graph
function create_topo_graph (topo_graph, gN, gC, resid_1, resid_2) {
  var  svgns = "http://www.w3.org/2000/svg";

  var center_x = 32
  var center_y = 28
  var r_1 = 12
  var r_2 = 24
  var angle_start = 45/180*3.14
  var angle_end = 170/180*3.14

  topo_graph.innerHTML = ""

  var graph = document.createElementNS(svgns, "svg")
  graph.setAttribute("width", "100%")
  graph.setAttribute("height", "100%")
  graph.setAttribute("viewBox", "0 0 64 64")
  
  var graph_html_list = []
  var curve_width = "2.5%"
  var edge_width = "5%"
  // add the closed loop
  var curve_color = "red"
  var edge_color = "white"
  var base_string = "<path d='M28 50 C28 48 28 46 24 44 C10 38 8 12 32 12 C56 12 54 38 40 44 C36 46 36 48 36 50' "
  graph_html_list.push(base_string+"stroke='"+edge_color+"' stroke-width='"+edge_width+"' fill='none'/>"+base_string+"stroke='"+curve_color+"' stroke-width='"+curve_width+"' fill='none'/>")

  var curve_color = "gray"
  var edge_color = "white"
  var base_string = "<path d='M22 52 C24 56 28 54 28 50' "
  graph_html_list.push(base_string+"stroke='"+edge_color+"' stroke-width='"+edge_width+"' fill='none'/>"+base_string+"stroke='"+curve_color+"' stroke-width='"+curve_width+"' fill='none'/>")
  var base_string = "<path d='M36 50 C36 54 40 56 42 52' "
  graph_html_list.push(base_string+"stroke='"+edge_color+"' stroke-width='"+edge_width+"' fill='none'/>"+base_string+"stroke='"+curve_color+"' stroke-width='"+curve_width+"' fill='none'/>")

  // add N-ter entanglements
  var end_x = 22;
  var end_y = 52;
  if (gN != 0) {
    var curve_color = "blue"
    var edge_color = "white"
    var delta_angle = (angle_end - angle_start) / Math.abs(gN);
    var i;
    for (i = 0; i < Math.abs(gN); i++) {
      var angle = angle_start + delta_angle * i;
      var x1 = center_x - r_1 * Math.sin(angle)
      var y1 = center_y + r_1 * Math.cos(angle)
      var x2 = center_x - r_2 * Math.sin(angle)
      var y2 = center_y + r_2 * Math.cos(angle)
      var base_string_1 = "<path d='M"+end_x.toString()+" "+end_y.toString()+" C"+(end_x-2).toString()+" "+(end_y-4).toString()+" "+(x1+2).toString()+" "+(y1+4).toString()+" "+x1.toString()+" "+y1.toString()+"' "
      var base_string_2 = "<path d='M"+x1.toString()+" "+y1.toString()+" C"+(x1-2).toString()+" "+(y1-4).toString()+" "+(x2+2).toString()+" "+(y2+4).toString()+" "+x2.toString()+" "+y2.toString()+"' "
      
      if (gN > 0) {
        graph_html_list.unshift(base_string_1+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string_1+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
        graph_html_list.push(base_string_2+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string_2+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
      }
      else {
        graph_html_list.push(base_string_1+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string_1+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
        graph_html_list.unshift(base_string_2+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string_2+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
      }

      end_x = x2
      end_y = y2
    }
  }
  var curve_color = "gray"
  var edge_color = "white"
  var x1 = 4
  var y1 = end_y
  var base_string = "<path d='M"+end_x.toString()+" "+end_y.toString()+" C"+(end_x-2).toString()+" "+(end_y-4).toString()+" "+(x1+4).toString()+" "+(y1+2).toString()+" "+x1.toString()+" "+y1.toString()+"' "
  graph_html_list.push(base_string+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
  graph_html_list.push("<text text-anchor='start' fill='black' font-size='30%' font-family='Arial' x=1 y="+(end_y+6).toString()+"> N-ter </text>")

  // add C-ter entanglements
  var end_x = 42;
  var end_y = 52;
  if (gC != 0) {
    var curve_color = "blue"
    var edge_color = "white"
    var delta_angle = (angle_end - angle_start) / Math.abs(gC);
    var i;
    for (i = 0; i < Math.abs(gC); i++) {
      var angle = angle_start + delta_angle * i;
      var x1 = center_x + r_1 * Math.sin(angle)
      var y1 = center_y + r_1 * Math.cos(angle)
      var x2 = center_x + r_2 * Math.sin(angle)
      var y2 = center_y + r_2 * Math.cos(angle)
      var base_string_1 = "<path d='M"+end_x.toString()+" "+end_y.toString()+" C"+(end_x+2).toString()+" "+(end_y-4).toString()+" "+(x1-2).toString()+" "+(y1+4).toString()+" "+x1.toString()+" "+y1.toString()+"' "
      var base_string_2 = "<path d='M"+x1.toString()+" "+y1.toString()+" C"+(x1+2).toString()+" "+(y1-4).toString()+" "+(x2-2).toString()+" "+(y2+4).toString()+" "+x2.toString()+" "+y2.toString()+"' "
      
      if (gC < 0) {
        graph_html_list.unshift(base_string_1+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string_1+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
        graph_html_list.push(base_string_2+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string_2+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
      }
      else {
        graph_html_list.push(base_string_1+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string_1+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
        graph_html_list.unshift(base_string_2+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string_2+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
      }

      end_x = x2
      end_y = y2
    }
  }
  var curve_color = "gray"
  var edge_color = "white"
  var x1 = 60
  var y1 = end_y
  var base_string = "<path d='M"+end_x.toString()+" "+end_y.toString()+" C"+(end_x+2).toString()+" "+(end_y-4).toString()+" "+(x1-4).toString()+" "+(y1+2).toString()+" "+x1.toString()+" "+y1.toString()+"' "
  graph_html_list.push(base_string+"stroke='"+edge_color+"' stroke-width='"+edge_width+"', fill='none'/>"+base_string+"stroke='"+curve_color+"' stroke-width='"+curve_width+"', fill='none'/>")
  graph_html_list.push("<text text-anchor='end' fill='black' font-size='30%' font-family='Arial' x=63 y="+(end_y+6).toString()+"> C-ter </text>")

  // add native contacts
  graph_html_list.push("<line x1=28 y1=50 x2=36 y2=50 stroke='orange' stroke-width='"+curve_width.toString()+"' stroke-dasharray='1% 0.5%' />")
  graph_html_list.push("<circle cx=28 cy=50 r=2 stroke='black' stroke-width='1%' fill='orange' />")
  graph_html_list.push("<circle cx=36 cy=50 r=2 stroke='black' stroke-width='1%' fill='orange' />")
  graph_html_list.push("<text text-anchor='end' fill='black' font-size='25%' font-family='Arial' x=30 y=60> "+resid_1.toString()+" </text>")
  graph_html_list.push("<text text-anchor='start' fill='black' font-size='25%' font-family='Arial' x=34 y=60> "+resid_2.toString()+" </text>")
  
  // add title
  var ggN = gN.toString()
  if (gN > 0) { ggN = "+"+ggN }
  var ggC = gC.toString()
  if (gC > 0) { ggC = "+"+ggC }
  graph_html_list.push("<text text-anchor='middle' fill='black' font-size='30%' font-family='Arial' x=32 y=6> "+
    "<tspan> <tspan font-style='italic'>g</tspan><tspan font-size='70%' baseline-shift='sub'>N</tspan> = "+ggN+"; <tspan font-style='italic'>g</tspan><tspan font-size='70%' baseline-shift='sub'>C</tspan> = "+ggC+" </tspan>"
    +" </text>")
  
  var graph_html = graph_html_list.join(" ")
  graph.innerHTML = graph_html

  topo_graph.appendChild(graph)
}

// load structure and show entanglements
var sele_dict = {}
var struct = []

function loadStructure (arg_str) {
  stage.removeAllComponents()

  var argvs = arg_str.split('&')
  var url_1 = "https://cdn.jsdelivr.net/gh/yuj179/topo_entanglements/structure_pdb/"+argvs[0]+".pdb" //url for state structure
  var s1 = "1-"+argvs[1] //select protein
  var s2 = argvs[2] //select circle
  var s3 = argvs[3] //select thread

  var words = s2.split(' or ')
  var s4 = ""
  for (i = 0; i < words.length; i++) {
    var w = words[i].split('(').slice(-1)[0].split(')')[0].split('-')
    s4 += w[0] + ' or ' + w[1] + ' or '
  }
  s4 = "(" + s4.substr(0, s4.length - 4) + ") and .CA" //select NC atoms

  var atom_pairs = [] //select NC atom pairs
  for (i = 0; i < words.length; i++) {
    var w = words[i].split('(').slice(-1)[0].split(')')[0].split('-')
    atom_pairs.push([w[0] + '.CA', w[1] + '.CA'])
  }

  var url_2 = "rcsb://"+argvs[4] //url for native protein pdb
  var chain_sel = ":"+argvs[5]
  var offset = Number(argvs[6])
  var s5 = chain_sel + " and " + (1+offset).toString()+"-"+(Number(argvs[1])+offset).toString() //select protein
  
  var s6 = chain_sel + " and " //select circle
  for (i = 0; i < words.length; i++) {
    var w = words[i].split('(').slice(-1)[0].split(')')[0].split('-')
    s6 += "("+(Number(w[0])+offset).toString()+"-"+(Number(w[1])+offset).toString()+") or "
  }
  s6 = s6.substr(0, s6.length - 4)

  var s7 = chain_sel + " and " //select thread
  words = s3.split(' or ')
  for (i = 0; i < words.length; i++) {
    var w = words[i].split('(').slice(-1)[0].split(')')[0].split('-')
    s7 += "("+(Number(w[0])+offset).toString()+"-"+(Number(w[1])+offset).toString()+") or "
  }
  s7 = s7.substr(0, s7.length - 4)

  var s8 = chain_sel + " and " + argvs[7] //select ligands
  
  var superpose_sel_1 = "(" + argvs[8]+") and .CA" //select for superpose of state structure
  var superpose_sel_2 = "" //select for superpose of state structure
  words = argvs[8].split(' or ')
  for (i = 0; i < words.length; i++) {
    var w = words[i].split('(').slice(-1)[0].split(')')[0].split('-')
    superpose_sel_2 += "("+(Number(w[0])+offset).toString()+"-"+(Number(w[1])+offset).toString()+") or "
  }
  superpose_sel_2 = "(" + superpose_sel_2.substr(0, superpose_sel_2.length - 4) + ") and .CA"
  //console.log(superpose_sel_1)

  var Q = argvs[9]
  var G = argvs[10]

  var if_ribosome = false
  if (Q == "-" && G == "-") {if_ribosome = true} 

  document.getElementById("info").innerHTML = "Loading structures... (if no response for a long time, please refresh this page)"
  
  // re-initialize the topology graph panels
  topo_graph_div.innerHTML = ""
  topo_graph_div.style.height = (div_length*atom_pairs.length).toString()+"vw"
  var i;
  var topo_graph_div_list = []
  for(i = 0; i < atom_pairs.length; i++) {
    var topo_sub_div = createElement("div", 
        { id: "topo_div_"+(i+1).toString() }, 
        { top: (i/atom_pairs.length*100)+"%", 
          left: "0%", 
          width: "100%",
          height: (1/atom_pairs.length*100)+"%",
          textAlign: "center",
          backgroundColor: "white",
          border: "2px solid black"})
    if (i != 0) {
      topo_sub_div.style.borderTop = "none"
    }
    topo_graph_div.appendChild(topo_sub_div)
    topo_graph_div_list.push(topo_sub_div)
  }
  
  // load structure and show entanglements
  Promise.all([
    stage.loadFile(url_1, {ext: 'pdb'}).then(function (o) {
      o.autoView()
      o.addRepresentation('cartoon', {
        sele: s1,
        color: "#E6E6FA",
        name: "protein_1"
      })
      sele_dict["protein_1"] = s1
      o.addRepresentation('cartoon', {
        sele: s2,
        color: "#FF0000",
        name: "circle_1"
      })
      sele_dict["circle_1"] = s2
      o.addRepresentation('cartoon', {
        sele: s3,
        color: "#0000FF",
        name: "thread_1"
      })
      sele_dict["thread_1"] = s3
      o.addRepresentation('spacefill', {
        sele: s4,
        color: "#FFA500",
        name: "NC_atoms"
      })
      o.addRepresentation('distance', {
        atomPair: atom_pairs,
        labelVisible: false,
        radius: 0.30,
        color: "#FFA500",
        name: "native contact"
      })
      return o
    }),

    stage.loadFile(url_2, {sele: chain_sel}).then(function (o) {
      o.autoView()
      o.addRepresentation('cartoon', {
        sele: s5,
        color: "#F0FFF0",
        name: "protein_2"
      })
      sele_dict["protein_2"] = s5
      o.addRepresentation('cartoon', {
        sele: s6,
        color: "#FFFF00",
        name: "circle_2"
      })
      sele_dict["circle_2"] = s6
      o.addRepresentation('cartoon', {
        sele: s7,
        color: "#32CD32",
        name: "thread_2"
      })
      sele_dict["thread_2"] = s7
      o.addRepresentation('spacefill', {
        sele: s8,
        scale: 0.8,
        name: "ligand"
      })
      sele_dict["ligand"] = s8
      return o
    })
  ]).then(function (ol) {
    ol[0].superpose(ol[1], false, superpose_sel_1, superpose_sel_2)
    ol[0].autoView()

    if (if_ribosome) {
        ol[0].addRepresentation("surface", {
          sele: "not ("+s1+" or :A or :P)",
          probeRadius: 6.0,
          scaleFactor: 0.35,
          surfaceType: "av",
          color: "white",
          opacity: 0.75,
          side: "front",
          opaqueBack: true,
          metalness: 0,
          roughness: 1,
          useWorker: false,
          name: "ribosome"
        })
        sele_dict["ribosome"] = "not ("+s1+" or :A or :P)"
        ol[0].addRepresentation("surface", {
          sele: "not ("+s1+") and (:A or :P)",
          probeRadius: 4.0,
          scaleFactor: 0.5,
          surfaceType: "av",
          color: "orange",
          metalness: 0,
          roughness: 1,
          useWorker: false,
          name: "tRNAs"
        })
        sele_dict["tRNAs"] = "not ("+s1+") and (:A or :P)"
      }
      else {
        sele_dict["ribosome"] = ""
        sele_dict["tRNAs"] = ""
      }

    struct = ol

    reset_view()

    document.getElementById("info").innerHTML = "Done loading structures. <i>Q</i><sub>act</sub> = "+Q+" and <i>G</i> = "+G

    // Get Ca trace
    var CA_cor_list = [];
    ol[0].structure.eachAtom(function (ap) {
      CA_cor_list.push([ap.x, ap.y, ap.z])
    }, new NGL.Selection(".CA"))
    // Get curve coordinates
    var R = []
    var dR = []
    var n_atom = CA_cor_list.length;
    var i;
    for(i = 0; i < n_atom-1; i++) {
      R.push([(CA_cor_list[i][0]+CA_cor_list[i+1][0])/2, (CA_cor_list[i][1]+CA_cor_list[i+1][1])/2, (CA_cor_list[i][2]+CA_cor_list[i+1][2])/2])
      dR.push([CA_cor_list[i+1][0]-CA_cor_list[i][0], CA_cor_list[i+1][1]-CA_cor_list[i][1], CA_cor_list[i+1][2]-CA_cor_list[i][2]])
    }
    // Gauss integral contents
    var M = [];
    var i;
    for(i = 0; i < n_atom-1; i++) {
      M.push([])
      var j;
      for(j = 0; j < n_atom-1; j++) {
        M[i].push(0)
      }
    }
    var i;
    for(i = 0; i < n_atom-2; i++) {
      var j;
      for(j = i+1; j < n_atom-1; j++) {
        var v1 = [R[i][0]-R[j][0], R[i][1]-R[j][1], R[i][2]-R[j][2]]
        var v1_n = (v1[0]**2+v1[1]**2+v1[2]**2)**(3/2)
        var v11 = [v1[0]/v1_n, v1[1]/v1_n, v1[2]/v1_n]
        var v2 = [dR[i][1]*dR[j][2]-dR[i][2]*dR[j][1], dR[i][2]*dR[j][0]-dR[i][0]*dR[j][2], dR[i][0]*dR[j][1]-dR[i][1]*dR[j][0]]
        M[i][j] = v11[0]*v2[0] + v11[1]*v2[1] + v11[2]*v2[2]
        M[j][i] = M[i][j]
      }
    }
    //Calculate gN and gC
    var i;
    for(i = 0; i < atom_pairs.length; i++) {
      var res_idx = []
      var sel = atom_pairs[i].join(" or ")
      ol[0].structure.eachAtom(function (ap) {
        res_idx.push(ap.resno)
      }, new NGL.Selection(sel))
      var g_list = []
      calc_g(g_list, M, res_idx[0], res_idx[1])
      create_topo_graph(topo_graph_div_list[i], g_list[0], g_list[1], res_idx[0], res_idx[1])
    }
  })
}

//////////////////////// Generate page contents //////////////////////////////
addElement(document.getElementById("cntrl_panel"), createElement("span", {
  innerText: "Select an example below:"
}, { top: "1%", 
     left: "1%", 
     "font-family": "Arial", 
     "font-size": "1.8vmin"}))

// Select syntax: 
// Github file name & chain length & circle_sel & thread_sel & Pdb ID & chain ID & offset & ligand_sel & superpose_sel & Q & G 
var StructSelect = createSelect([
  [ "CAT-III_state_P12&213&167-183&25-44&3cla&A&5&CLM&1-213&0.87&0.13", "CAT-III State P12" ],
  [ "CAT-III_state_P10&213&85-101&122-141&3cla&A&5&CLM&1-213&0.64&0.05", "CAT-III State P10" ],
  [ "CAT-III_state_P13&213&71-163&167-186&3cla&A&5&CLM&1-213&0.85&0.09", "CAT-III State P13" ],
  [ "DDLB_state_P5&306&119-177&185-204&4c5c&A&0&(310-312) or (1313-1314)&1-306&0.66&0.06", "DDLB State P5" ],
  [ "DDLB_state_P10&306&145-174&121-140&4c5c&A&0&(310-312) or (1313-1314)&1-306&0.85&0.04", "DDLB State P10" ],
  [ "DDLB_state_C8&306&119-177&179-198&4c5c&A&0&(310-312) or (1313-1314)&119-177&-&-", "DDLB State C8" ],
  [ "DDLB_state_C9&306&145-174&121-140&4c5c&A&0&(310-312) or (1313-1314)&145-174&-&-", "DDLB State C9" ],
], {
  onchange: function (e) {
    loadStructure(e.target.value)
  }
}, { top: "4%", 
     left: "1%",
     "font-family": "Arial", 
     "font-size": "1.8vmin"})
addElement(document.getElementById("cntrl_panel"), StructSelect)

addElement(document.getElementById("cntrl_panel"), createElement("span", {
  innerText: "Check/uncheck the following \noptions to show/hide:"
}, { top: "7%", 
     left: "1%",
     "font-family": "Arial", 
     "font-size": "1.8vmin"}))

var ligandCheckbox = createElement("input", {
  type: "checkbox",
  checked: true,
  onchange: function (e) {
    stage.getRepresentationsByName("ligand")
      .setVisibility(e.target.checked)
  }
}, { top: "12%", 
     left: "1%", 
     width: "1.6vmin",
     height: "1.6vmin",
     margin: "0.2vmin" })
addElement(document.getElementById("cntrl_panel"), ligandCheckbox)
addElement(document.getElementById("cntrl_panel"), createElement("span", {
  innerText: "ligand"
}, { top: "12%", 
     left: "5vmin",
     "font-family": "Arial", 
     "font-size": "1.8vmin" }))

var nativeCheckbox = createElement("input", {
  type: "checkbox",
  checked: true,
  onchange: function (e) {
    stage.getRepresentationsByName("protein_2")
      .setVisibility(e.target.checked)
    stage.getRepresentationsByName("circle_2")
      .setVisibility(e.target.checked)
    stage.getRepresentationsByName("thread_2")
      .setVisibility(e.target.checked)
  }
}, { top: "14.5%", 
     left: "1%", 
     width: "1.6vmin",
     height: "1.6vmin",
     margin: "0.2vmin" })
addElement(document.getElementById("cntrl_panel"), nativeCheckbox)
addElement(document.getElementById("cntrl_panel"), createElement("span", {
  innerText: "native structure"
}, { top: "14.5%", 
     left: "5vmin",
     "font-family": "Arial", 
     "font-size": "1.8vmin" }))

var stateCheckbox = createElement("input", {
  type: "checkbox",
  checked: true,
  onchange: function (e) {
    stage.getRepresentationsByName("protein_1")
      .setVisibility(e.target.checked)
    stage.getRepresentationsByName("ribosome")
      .setVisibility(stateCheckbox.checked)
    stage.getRepresentationsByName("tRNAs")
      .setVisibility(stateCheckbox.checked)
  }
}, { top: "17%", 
     left: "1%", 
     width: "1.6vmin",
     height: "1.6vmin",
     margin: "0.2vmin" })
addElement(document.getElementById("cntrl_panel"), stateCheckbox)
addElement(document.getElementById("cntrl_panel"), createElement("span", {
  innerText: "misfolding structure"
}, { top: "17%", 
     left: "5vmin",
     "font-family": "Arial", 
     "font-size": "1.8vmin" }))

var entanglementCheckbox = createElement("input", {
  type: "checkbox",
  checked: true,
  onchange: function (e) {
    stage.getRepresentationsByName("circle_1")
      .setVisibility(e.target.checked)
    stage.getRepresentationsByName("thread_1")
      .setVisibility(e.target.checked)
    stage.getRepresentationsByName("NC_atoms")
      .setVisibility(e.target.checked)
    stage.getRepresentationsByName("native contact")
      .setVisibility(e.target.checked)
  }
}, { top: "19.5%", 
     left: "1%", 
     width: "1.6vmin",
     height: "1.6vmin",
     margin: "0.2vmin" })
addElement(document.getElementById("cntrl_panel"), entanglementCheckbox)
addElement(document.getElementById("cntrl_panel"), createElement("span", {
  innerText: "entanglements"
}, { top: "19.5%", 
     left: "5vmin",
     "font-family": "Arial", 
     "font-size": "1.8vmin" }))

var centerButton = createElement("input", {
  type: "button",
  value: "center",
  onclick: function () {
    reset_view()
  }
}, { top: "22.5%", 
     left: "1%",
     "font-family": "Arial", 
     "font-size": "1.8vmin" })
addElement(document.getElementById("cntrl_panel"), centerButton)

addElement(document.getElementById("cntrl_panel"), createElement("span", {
  innerText: "Lasso-like entanglement topology:"
}, { top: "27%", 
     left: "1%", 
     "font-family": "Arial", 
     "font-size": "1.8vmin"}))

// initialize the topology graph panels
var div_length = 15
var topo_graph_div = createElement("div", 
  { id: "topo_div" }, 
  { top: "30%", 
    left: "1%", 
    width: div_length.toString()+"%",
    height: div_length.toString()+"vw",
    textAlign: "center",
    backgroundColor: "white"})
addElement(document.getElementById("cntrl_panel"), topo_graph_div)

loadStructure("CAT-III_state_P12&213&167-183&25-44&3cla&A&5&CLM&1-213&0.87&0.13")
