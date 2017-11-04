using UnityEngine;
using System.Collections;
using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

public class findSurface : MonoBehaviour {

	// Use this for initialization

    //Mesh m = new Mesh();

    GameObject tet;
    MeshFilter meshFilter;
    public static Matrix<double> node;
    public static Matrix<double> element;


      

	void Start () {

        

        node = readNodes();
        element = readElements();

        int numElements = element.RowCount;
        int numNodes = node.RowCount;

        Debug.Log("Number of Elements: " + numElements);
        Debug.Log("Number of Nodes: " + numNodes);


        

        //4 choose 3 to create 4 faces for each of the elements
        Matrix<double> faces = buildFaces(element, numElements);

        //sort each row in ascending order horizontally

        faces = sortRows(faces, numElements);

       

        Matrix<double> surface = determineSurfaceFaces(faces, numElements);

        

        int numFaces = surface.RowCount;

        //Vector<double> surfaceTets = determineSurfaceTets(faces,surface, numFaces);



        //create mesh

       

        Vector3[] vertices = new Vector3[numNodes];

        for (int i = 0; i < numNodes; i++)  //fill in the vertices matrix
        {
            vertices[i] = new Vector3((float)(node[i, 0]), (float)(node[i, 1]), (float)(node[i, 2]));
        }


        int[] triangles = new int[numFaces * 3];

        int z = 0;

        for (int j = 0; j < numFaces; j++)
        {
            for (int q = 0; q < 3; q++)
            {
                triangles[z] = (int)(surface[j, q]) - 1;
                z++;
            }
        }

 

        triangles = flipFaces(surface, numFaces, triangles);


        


        Mesh m = new Mesh();


        m.name = "ScriptedMesh";
        

        m.vertices = vertices;
        m.triangles = triangles;
        
        tet = new GameObject("Tetahedral Sphere Surface");
        meshFilter = (MeshFilter)tet.AddComponent(typeof(MeshFilter));

        Texture2D tex = new Texture2D(1, 1);
        MeshRenderer renderer = tet.AddComponent(typeof(MeshRenderer)) as MeshRenderer;
        ;
        meshFilter.mesh = m;
        renderer.material.shader = Shader.Find("Standard");
        
        tex.SetPixel(0, 0, Color.green);
        tex.Apply();
        //renderer.material.mainTexture = tex;
        //renderer.material.color = Color.green;

        //renderer.material = thisColor;
        // Assigns a material named "Assets/Resources/DEV_Orange" to the object.
        Material newMat = Resources.Load("thisColor", typeof(Material)) as Material;  
        Texture newText = Resources.Load("vivemap", typeof(Texture)) as Texture;  
        renderer.material = newMat;
        renderer.material.mainTexture = newText;

    
        

        
    }

    // Update is called once per frame
    void Update()
    {

        //int numElements = element.RowCount;
        //int numNodes = node.RowCount;

        if (Input.GetKeyDown("d"))
        {
            /*

            


            //element.RemoveRow(surfaceTets[i]);
            numElements--;


            //4 choose 3 to create 4 faces for each of the elements
            Matrix<double> faces = buildFaces(element, numElements);

            //sort each row in ascending order horizontally

            faces = sortRows(faces, numElements);



            Matrix<double> surface = determineSurfaceFaces(faces, numElements);

            int numFaces = surface.RowCount;

            Vector<double> surfaceTets = determineSurfaceTets(faces, surface, numFaces);
            //create mesh




            Vector3[] vert = new Vector3[numNodes];

            for (int i = 0; i < numNodes; i++)  //fill in the vertices matrix
            {
                vert[i] = new Vector3((float)(node[i, 0]), (float)(node[i, 1]), (float)(node[i, 2]));
            }


            int[] tri = new int[surface.RowCount * 3];

            int z = 0;

            for (int j = 0; j < surface.RowCount; j++)
            {
                for (int q = 0; q < 3; q++)
                {
                    tri[z] = (int)(surface[j, q]) - 1;
                    z++;
                }
            }

            Vector3[] vertices = new Vector3[numNodes];

            for (int i = 0; i < numNodes; i++)  //fill in the vertices matrix
            {
                vertices[i] = new Vector3((float)(node[i, 0]), (float)(node[i, 1]), (float)(node[i, 2]));
            }

            //tri = flipFaces(surface, surface.RowCount, tri);

            Mesh m = new Mesh();

            m.vertices = vert;
            m.triangles = tri;

            //m.RecalculateNormals();
            meshFilter.mesh = m;


            Debug.Log(element);
            Debug.Log("Row Count: " + element.RowCount);*/
        }

    }

    public static Matrix<double> readNodes()
    {
        
        string line;
        List<double> nodes = new List<double>();

     

        System.IO.StreamReader nodes_file = new System.IO.StreamReader("./assets/sphere_nodes.txt");
        while ((line = nodes_file.ReadLine()) != null)
        {

            char[] delimiterChars = { ' ',','};

            //Debug.Log(line);
            string[] node_chars = line.Split(delimiterChars);
            int n;
            foreach (string s in node_chars)
            {
               
                    nodes.Add(Convert.ToDouble(s));
            
            }

        }

        nodes_file.Close();
        //Debug.Log(nodes.Count/3);
        Matrix<double> node = Matrix<double>.Build.Dense(nodes.Count/3, 3);

        int counter = 0;

        for (int i = 0; i < (nodes.Count / 3); i++)
        {
            for (int j = 0; j < 3; j++)
            {
        

                    node[i, j] = nodes[counter];
                counter++;
            }
        }

        return node;

    }

    public static Matrix<double> readElements()
    {

        string line;
        List<double> elements = new List<double>();

        System.IO.StreamReader element_file = new System.IO.StreamReader(@"./assets/sphere_elements.txt");
        while ((line = element_file.ReadLine()) != null)
        {

            char[] delimiterChars = { ' ', ',','N' };

            string[] element_chars = line.Split(delimiterChars);

            foreach (string s in element_chars)
            {
                
                    elements.Add(Convert.ToDouble(s));
                
            }

        }

        element_file.Close();

        Matrix<double> element = Matrix<double>.Build.Dense(elements.Count / 4, 4);

        int counter = 0;

        for (int i = 0; i < (elements.Count / 4); i++)
        {
            for (int j = 0; j < 4; j++)
            {

                

                element[i, j] = elements[counter];
                counter++;
            }
        }

        return element;

    }


    public double TetrahedronElementVolume(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4)
    {
        Matrix<double> xyz;
        double V;
        xyz = DenseMatrix.OfArray(new double[,] {
			{1.0, x1, y1, z1},
			{1, x2, y2, z2},
			{1.0, x3, y3, z3},
			{1.0, x4, y4, z4}});
        V = xyz.Determinant() / 6;
        return V;
    }

    
    public Matrix<double> buildFaces(Matrix<double> element, int numElements)
    {
        Matrix<double> faces = Matrix<double>.Build.Dense(numElements*4, 3);
        
        int k=0;

        for (int i = 0; i < numElements; i++)
        {

            faces[k, 0] = element[i, 0];
            faces[k, 1] = element[i, 1];
            faces[k, 2] = element[i, 2];
            k++;
        }

        for (int i = 0; i < numElements; i++)
        {

            faces[k, 0] = element[i, 0];
            faces[k, 1] = element[i, 1];
            faces[k, 2] = element[i, 3];
            k++;
        }

        for (int i = 0; i < numElements; i++)
        {

            faces[k, 0] = element[i, 0];
            faces[k, 1] = element[i, 2];
            faces[k, 2] = element[i, 3];
            k++;
        }

        for (int i = 0; i < numElements; i++)
        {

            faces[k, 0] = element[i, 1];
            faces[k, 1] = element[i, 2];
            faces[k, 2] = element[i, 3];
            k++;
        }



        return faces;
	}

    public Matrix<double> sortRows(Matrix<double> faces, int numElements)
    {

        double temp = 0.0;

        for (int i = 0; i < numElements*4; i++)
        {

            if (faces[i, 0] > faces[i, 2])
            {
                temp = faces[i, 0];
                faces[i, 0] = faces[i, 2];
                faces[i, 2] = temp;
            }

            if (faces[i, 0] > faces[i, 1])
            {
                temp = faces[i, 0];
                faces[i, 0] = faces[i, 1];
                faces[i, 1] = temp;
            }

            if (faces[i, 1] > faces[i, 2])
            {
                temp = faces[i, 1];
                faces[i, 1] = faces[i, 2];
                faces[i, 2] = temp;
            }

        }

        return faces;

    }

    public Matrix<double> determineSurfaceFaces(Matrix<double> faces, int numElements)
    {


        HashSet<Vector<double>> surfaceSet = new HashSet<Vector<double>>();
        
        for (int i = 0; i < numElements*4; i++)
        {
            if (surfaceSet.Contains(faces.Row(i)))
                surfaceSet.Remove(faces.Row(i));
            else
                surfaceSet.Add(faces.Row(i));
        }


       Matrix<double> surface = Matrix<double>.Build.Dense(surfaceSet.Count, 3);

       int row = 0;
       foreach (Vector<double> temp in surfaceSet)
       {
           surface.SetRow(row, temp);
           row++;
       }

           return surface;


    }




    public int[] flipFaces(Matrix<double> f, int numFaces, int[] triangles)
    {

        System.Diagnostics.Stopwatch stopWatch = new System.Diagnostics.Stopwatch();

        int numEdges = numFaces * 3;


        
        


        //create 3D array with edge indices
        double[, ,] edges = create3Dedges(f, numFaces);

     

        //mask which indicates which edges have indicies in ascending order
        Matrix<double> edgeAscOrderMask = ascending(edges, numFaces);

      

        //flatten edges into (numFaces*3)x2 matrix
        Matrix<double> edgesFlat = flatten(edges, numEdges, numFaces);

        stopWatch.Reset();
        stopWatch.Start();

        //determine edge group numbers
        Vector<double> edgeGrpNos = Vector<double>.Build.Dense(numFaces * 3);
        int grpNos = groupEdges(edgeGrpNos, edgesFlat, numFaces);

        stopWatch.Stop();
        //Debug.Log("Function 4: " + stopWatch.Elapsed.Milliseconds);
        stopWatch.Reset();
        stopWatch.Start();

        //Determine which edges are nicely partnered (asc/desc), and the face
        //number that they link to
        Vector<double> edgeGrpIsUnified = Vector<double>.Build.Dense(numEdges);
        Vector<double> edgePartnerFaceNo = Vector<double>.Build.Dense(numEdges);
        unifiedPartners(edgeGrpIsUnified, edgePartnerFaceNo, numEdges, numFaces, grpNos, edgeGrpNos, edgeAscOrderMask);

        stopWatch.Stop();
        //Debug.Log("Function 5: " + stopWatch.Elapsed.Milliseconds);
        stopWatch.Reset();
        stopWatch.Start();


        //Collect connected faces/edges
        //March from the first face to each of its nicely (asc/desc) connected
        //neighbour faces. Label each connected "set" of faces.
        Vector<double> faceSets = Vector<double>.Build.Dense(numFaces);
        int currentSet = groupFaces(faceSets, edgeGrpIsUnified, edgePartnerFaceNo, numEdges, numFaces);

        stopWatch.Stop();
        //Debug.Log("Function 6: " + stopWatch.Elapsed.Milliseconds);
        stopWatch.Reset();
        stopWatch.Start();


        //Determine which sets need to be flipped
        Vector<double> setsTouched = Vector<double>.Build.Dense(currentSet);
        Vector<double> setsToFlip = Vector<double>.Build.Dense(currentSet);
        whichMustFlip(setsTouched, setsToFlip, faceSets, edgePartnerFaceNo, numFaces, currentSet);


        stopWatch.Stop();
        //Debug.Log("Function 7: " + stopWatch.Elapsed.Milliseconds);
        stopWatch.Reset();
        stopWatch.Start();

        //Perform the actual flipping of all sets that require it
        //Debug.Log(f);

        int[] newTriangles = flipping(setsToFlip, faceSets, triangles, numEdges, numFaces);


        stopWatch.Stop();
        //Debug.Log("Function 8: " + stopWatch.Elapsed.Milliseconds);
        stopWatch.Reset();
        stopWatch.Start();

        return newTriangles;


    }


    public double[, ,] create3Dedges(Matrix<double> f, int numFaces)
    {

        double[, ,] edges = new double[numFaces, 3, 2];//create 3D array with edge indices

        for (int i = 0; i < 2; i++)
        {//3rd dimension
            if (i == 0)
            {
                for (int j = 0; j < numFaces; j++)//rows
                {
                    for (int c = 0; c < 3; c++)//columns
                    {
                        edges[j, c, i] = f[j, c];
                    }
                }
            }
            else if (i == 1)
            {
                for (int j = 0; j < numFaces; j++)//rows
                {
                    for (int c = 0; c < 3; c++)//columns
                    {
                        if (c == 2)
                        {
                            edges[j, c, i] = f[j, c - 2];
                        }
                        else
                        {
                            edges[j, c, i] = f[j, c + 1];
                        }
                    }
                }
            }
        }//end of 3rd dimension

        return edges;

    }

    public Matrix<double> ascending(double[, ,] edges, int numFaces)
    {
        Matrix<double> edgeAscOrderMask = Matrix<double>.Build.Dense(numFaces, 3);

        for (int j = 0; j < numFaces; j++)//rows
        {
            for (int c = 0; c < 3; c++)//columns
            {
                if (edges[j, c, 0] < edges[j, c, 1])
                {
                    edgeAscOrderMask[j, c] = 1;

                }
                else if (edges[j, c, 0] == edges[j, c, 1])
                {

                }
                else
                {
                    edgeAscOrderMask[j, c] = 0;

                }

            }
        }

        return edgeAscOrderMask;
    }

    public Matrix<double> flatten(double[, ,] edges, int numEdges, int numFaces)
    {

        Matrix<double> edgesFlat = Matrix<double>.Build.Dense(numEdges, 2);

        for (int c = 0; c < 3; c++)//columns
        {
            for (int j = 0; j < numFaces; j++)
            {
                edgesFlat[c * numFaces + j, 0] = edges[j, c, 0];
                edgesFlat[c * numFaces + j, 1] = edges[j, c, 1];
            }

        }

        //swap the vertices so the smaller index for each edge is in first column
        for (int j = 0; j < numEdges; j++)
        {
            if (edgesFlat[j, 0] > edgesFlat[j, 1])
            {
                double holder = edgesFlat[j, 0];
                edgesFlat[j, 0] = edgesFlat[j, 1];
                edgesFlat[j, 1] = holder;
            }
        }
        return edgesFlat;

    }

    public int groupEdges(Vector<double> edgeGrpNos, Matrix<double> edgesFlat, int numFaces)
    {
        List<Vector<double>> edgeGrp = new List<Vector<double>>();
        int grpNos = 0;

        for (int j = 0; j < numFaces * 3; j++)
        {
            Vector<double> edge = edgesFlat.Row(j);
            if (!edgeGrp.Contains(edge))
            {
                edgeGrp.Add(edge);
                edgeGrpNos[j] = grpNos;
                grpNos++;

            }
            else
            {
                edgeGrpNos[j] = edgeGrp.IndexOf(edge);

            }
        }

        return grpNos;
    }

    public void unifiedPartners(Vector<double> edgeGrpIsUnified, Vector<double> edgePartnerFaceNo, int numEdges, int numFaces, int grpNos, Vector<double> edgeGrpNos, Matrix<double> edgeAscOrderMask)
    {

        
        Vector<double> mask = Vector<double>.Build.Dense(numEdges);
        Vector<double> mask3d = Vector<double>.Build.Dense(numEdges);
        Vector<double> inds = Vector<double>.Build.Dense(2);
        //Debug.Log(grpNos);

     


        for (int i = 0; i < grpNos; i++)
        {
            //create mask selecting row which contains index of current i
            mask.Clear();
            mask3d.Clear();

            for (int j = 0; j < 3; j++)
            {
                for (int w = 0; w < numFaces; w++)
                {
                    if (edgeGrpNos[numFaces * j + w] == i)
                    {
                        mask[numFaces * j + w] = 1;
                        mask3d[numFaces * j + w] = edgeAscOrderMask[w, j];
                    }
                }
            }

            //set new matrix edgeGrpIsUnified to 1 at mask if there are 2
            //occurances and if values at mask in edgeAscOrderMask add to 1
            if (mask.Sum() == 2.0 && mask3d.Sum() == 1.0)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int w = 0; w < numFaces; w++)
                    {
                        if (edgeGrpNos[numFaces * j + w] == i)
                            edgeGrpIsUnified[numFaces * j + w] = 1;
                    }
                }
            }

            

            //declare inds
            int firstfound = 0;
            for (int p = 0; p < numEdges; p++)
            {
                if (mask[p] == 1)
                {
                    if (firstfound == 0)
                    {
                        inds[0] = p;
                        firstfound = 1;
                    }
                    else
                    {
                        inds[1] = p;
                    }
                }
            }

            //Return which face each edge is part of
            for (int p = 0; p < numEdges; p++)
            {
                if (p == inds[0])
                {
                    edgePartnerFaceNo[p] = inds[1] % numFaces;
                }
                else if (p == inds[1])
                {
                    edgePartnerFaceNo[p] = inds[0] % numFaces;
                }
            }
        }
    }

    public int groupFaces(Vector<double> faceSets, Vector<double> edgeGrpIsUnified, Vector<double> edgePartnerFaceNo, int numEdges, int numFaces)
    {

        int currentSet = 0;
        Vector<double> facesLocked = Vector<double>.Build.Dense(numFaces);
        Vector<double> edgesConsidered = Vector<double>.Build.Dense(numEdges);
        List<int> currFaces = new List<int>();


        //decalare unit spaced
        Vector<double> unitSpaced = Vector<double>.Build.Dense(numEdges / numFaces);

        int k = 0;
        for (int i = 0; i < numEdges; i++)
        {
            if (i % numFaces == 0)
            {
                unitSpaced[k] = i;
                k++;
            }
        }

        // Debug.Log(unitSpaced);

        //decalare unit spaced


        while (facesLocked.Sum() != numFaces)
        {
            //If we're not connected to anything, we must start a new set
            if (currFaces.Count() == 0)
            {
                //find first non-zero element
                for (int i = 0; i < numFaces; i++)
                {
                    if (facesLocked[i] == 0)
                    {
                        currFaces.Add(i);
                        break;
                    }
                }

                currentSet++;
            }

            //lock the current face at indices of currFaces

            for (int i = 0; i < currFaces.Count(); i++)
            {
                for (int j = 0; j < numFaces; j++)
                {
                    if (j == currFaces[i])
                    {
                        facesLocked[j] = 1;
                        faceSets[j] = currentSet;
                    }
                }
            }




            //Grab the edges of the current faces
            Vector<double> facesVector = Vector<double>.Build.Dense(currFaces.Count()); //currFaces as a vector

            Vector<double> currEdgeInds = Vector<double>.Build.Dense(currFaces.Count() * unitSpaced.Count());

            //convert the list into a matrix for manipulation
            for (int i = 0; i < currFaces.Count(); i++)
            {
                facesVector[i] = currFaces[i];


            }

            k = 0;
            for (int i = 0; i < facesVector.Count(); i++)
            {
                for (int j = 0; j < unitSpaced.Count(); j++)
                {
                    currEdgeInds[k] = facesVector[i] + unitSpaced[j];
                    k++;
                }
            }




            //Find which edges are nicely connected to unvisited faces
            List<double> linked = new List<double>();

            for (int i = 0; i < currEdgeInds.Count(); i++)
            {
                if (edgeGrpIsUnified[(int)(currEdgeInds[i])] == 1 && edgesConsidered[(int)(currEdgeInds[i])] == 0)
                {

                    linked.Add(edgePartnerFaceNo[(int)(currEdgeInds[i])]);
                }

                //update the visited edges
                for (int j = 0; j < numEdges; j++)
                {
                    if (j == currEdgeInds[i])
                    {
                        edgesConsidered[j] = 1;
                    }
                }
            }

            //Determine the new faces we would reach if we stepped via nice edges
            if (linked.Count() > 0)
            {
                Vector<double> linkedFaces = Vector<double>.Build.Dense(linked.Count());

                for (int i = 0; i < linked.Count(); i++)
                {
                    linkedFaces[i] = linked[i];
                }


                currFaces.Clear();

                for (int i = 0; i < linkedFaces.Count(); i++)
                {

                    if (facesLocked[(int)(linkedFaces[i])] == 0)
                    {
                        currFaces.Add((int)(linkedFaces[i]));
                    }

                }
            }
            else
            {
                currFaces.Clear();
            }
        }
        return currentSet;
    }

    public void whichMustFlip(Vector<double> setsTouched, Vector<double> setsToFlip, Vector<double> faceSets, Vector<double> edgePartnerFaceNo, int numFaces, int currentSet)
    {

        List<double> currentSets = new List<double>();
        bool flipTheNextSet = true;

        while (setsTouched.Sum() != currentSet)
        {

            //If no current sets, pick the first one available. Any (next) sets
            //found touching it will need to be flipped

            if (currentSets.Count() == 0)
            {
                for (int i = 0; i < setsTouched.Count(); i++)
                {
                    if (setsTouched[i] == 0)
                    {
                        currentSets.Add(i + 1);   //add the set number if the set has not been touched yet (setNo VS. indexNo)
                        break;
                    }
                }
                flipTheNextSet = true;
            }



            //We've now touched the current sets. Find these sets' faces.
            for (int i = 0; i < setsTouched.Count(); i++)
            {
                for (int j = 0; j < currentSets.Count(); j++)
                {
                    if (i + 1 == currentSets[j]) //update setsTouched at index value in currentSet
                        setsTouched[i] = 1;
                }
            }

            //Find all the faces that share these edges

            List<double> setsFaceInds = new List<double>();

            List<double> unqFaceNosSharingBorder = new List<double>();


            for (int i = 0; i < faceSets.Count(); i++)
            {
                for (int j = 0; j < currentSets.Count(); j++)
                {
                    if (faceSets[i] == currentSets[j])
                    {
                        setsFaceInds.Add(i + 1);  //add the sets that are present in current set
                    }
                }
            }




            for (int i = 0; i < edgePartnerFaceNo.Count(); i++)
            {
                for (int j = 0; j < setsFaceInds.Count(); j++)
                {
                    if (edgePartnerFaceNo[i] == (setsFaceInds[j] - 1))
                    {
                        if (!unqFaceNosSharingBorder.Contains(i % numFaces))
                            unqFaceNosSharingBorder.Add(i % numFaces);

                    }
                }
            }






            List<double> setNosToFlip = new List<double>();

            if (unqFaceNosSharingBorder.Count() > 0)
            {



                //Avoid flipping faces on any sets already touched

                //We have a set of unqFaceNosSharing Border
                //Find out which set they are in and if it has been touched
                //if it has not, add it to faceNosToFlip




                for (int i = 0; i < unqFaceNosSharingBorder.Count(); i++)
                {
                    for (int j = 0; j < faceSets.Count(); j++)
                    {
                        if (unqFaceNosSharingBorder[i] == j)
                        {

                            if (setsTouched[(int)(faceSets[j]) - 1] == 0)
                            {
                                setNosToFlip.Add(faceSets[j]);
                            }

                        }
                    }
                }


            }
            else
            {
                setNosToFlip.Clear();
            }

            //Flip those sets if we should. The first (root) set WON'T have been
            //flipped. Any touching it WILL get flipped. Any touching *those* will
            //already be in the same direction as the root set, so they WON'T be
            //flipped. The next WILL, WON'T, WILL, WON'T, etc.


            if (flipTheNextSet)
            {
                for (int i = 0; i < setsToFlip.Count(); i++)
                {
                    for (int j = 0; j < setNosToFlip.Count(); j++)
                    {
                        if ((i + 1) == setNosToFlip[j]) //if setNoToFlip is at index location of sets to flip, make it true;
                        {
                            setsToFlip[i] = 1;
                        }

                    }
                }


            }
            flipTheNextSet = !flipTheNextSet;
            //Let the loop continue using the sets we just flipped as a source.

            currentSets.Clear();
            for (int i = 0; i < setNosToFlip.Count(); i++)
            {
                currentSets.Add(setNosToFlip[i]);
            }

            //check if there is a set that hasn't been touched yet




        }
    }

    public int[] flipping(Vector<double> setsToFlip, Vector<double> faceSets, int[] triangles, int numEdges, int numFaces)
    {


        int k = 0;

        Matrix<double> finalFaces = Matrix<double>.Build.Dense(numFaces, 3);

        for (int i = 0; i < numFaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                finalFaces[i, j] = triangles[k];
                k++;
            }
        }




        double placeholder;

        for (int i = 0; i < faceSets.Count(); i++)
        {
            for (int j = 0; j < setsToFlip.Count(); j++)
            {
                if (setsToFlip[j] == 1 && faceSets[i] == (j + 1)) //if the set must be flipped and the face corresponding to that set shows up, flip it
                {
                    placeholder = finalFaces[i, 1];
                    finalFaces[i, 1] = finalFaces[i, 0];
                    finalFaces[i, 0] = placeholder;
                }
            }
        }
       // Debug.Log(finalFaces);
        int numNodes = node.RowCount;
        
        Vector3[] v = new Vector3[numNodes];

        Vector<double> column_sums = node.ColumnSums();

        double x_mean = column_sums[0]/numNodes;
        double y_mean = column_sums[1]/numNodes;
        double z_mean = column_sums[2]/numNodes;
        

        for (int i = 0; i < numNodes; i++)  //fill in the vertices matrix
        {
            v[i] = new Vector3((float)(node[i, 0]-x_mean), (float)(node[i, 1]-y_mean), (float)(node[i, 2]-z_mean));
        }

        Matrix<double> tetra = Matrix<double>.Build.Dense(3, 3);
        Vector<double> tetvols = Vector<double>.Build.Dense((int)(element.RowCount));
        
        
        //perform outward alignment
        for(int i=0 ; i<numFaces ; i++){

           tetra[0,0] = v[(int)(finalFaces[i, 0])].x;
           tetra[0,1] = v[(int)(finalFaces[i, 0])].y;
           tetra[0,2] = v[(int)(finalFaces[i, 0])].z;

           tetra[1,0] = v[(int)(finalFaces[i, 1])].x;
           tetra[1,1] = v[(int)(finalFaces[i, 1])].y;
           tetra[1,2] = v[(int)(finalFaces[i, 1])].z;

           tetra[2,0] = v[(int)(finalFaces[i, 2])].x;
           tetra[2,1] = v[(int)(finalFaces[i, 2])].y;
           tetra[2,2] = v[(int)(finalFaces[i, 2])].z;

           tetvols[i] = tetra.Determinant() / 6;

        }

        double sum = tetvols.Sum();

       
        if (sum < 0)
        {
            for (int i = 0; i < numFaces; i++)
            {
                
                    placeholder = finalFaces[i, 1];
                    finalFaces[i, 1] = finalFaces[i, 0];
                    finalFaces[i, 0] = placeholder;
                
            }

        }

        
       


        //assign triangles array to new faces and assign to mesh

        int[] newTriangles = new int[numEdges];

        int z = 0;

        for (int j = 0; j < numFaces; j++)
        {
            for (int q = 0; q < 3; q++)
            {
                newTriangles[z] = (int)(finalFaces[j, q]);
                z++;
            }
        }
        return newTriangles;
    }

    public Matrix<double> cleanUpFaces(Vector3[] vertices, int[] triangles)
    {

        int numVertices = vertices.Length;
        int numFaces = triangles.Length / 3;



        Matrix<double> faces_orig = Matrix<double>.Build.Dense(numFaces, 3); //build matrix of faces
        Matrix<double> vertices_orig = Matrix<double>.Build.Dense(numVertices, 3); //build matrix of faces

        int k = 0;                            //fill the faces matrix
        for (int i = 0; i < numFaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                faces_orig[i, j] = triangles[k];
                k++;
            }

        }

        for (int i = 0; i < numVertices; i++)  //fill in the vertices matrix
        {
            vertices_orig[i, 0] = vertices[i].x;
            vertices_orig[i, 1] = vertices[i].y;
            vertices_orig[i, 2] = vertices[i].z;
        }



        /* CLEAN UP REPEATED VERTICES LABELS IN BOTH VERTICES AND FACES MATRICES 
         * 
         * 1.Determine indices of vertices matrix which contain duplicated vertices. Delete all but the first occurances.
         * 2.Replace those values of indicies with first index in faces matrix
           3.Create new faces and vertices matrices.
         */

        //create copy faces and vertices array
        Matrix<double> f = faces_orig;
        Matrix<double> v = vertices_orig;

        HashSet<double> currRepeats = new HashSet<double>();
        HashSet<int> totalRepeats = new HashSet<int>();


        //Debug.Log(v);
        //Debug.Log(f);
        //DelimitedWriter.Write("face_before.csv", f, " ");
        //DelimitedWriter.Write("vertex_before.csv", v, " ");

        for (int i = 0; i < numVertices; i++)
        {
            //store value from first row
            Vector<double> original = v.Row(i);
            //loop through the rest of the rows
            for (int j = i + 1; j < numVertices; j++)
            {
                //if the current row matches, save the index
                Vector<double> currRow = v.Row(j);

                if (currRow.Equals(original))
                {
                    currRepeats.Add(j);
                    totalRepeats.Add(j);
                }
            }
            //once finished looping, replace all the values in the faces array
            //move to next line and repeat

            for (int q = 0; q < numFaces; q++)
            {
                for (int w = 0; w < 3; w++)
                {
                    if (currRepeats.Contains(f[q, w]))
                        f[q, w] = i;
                }
            }

            currRepeats.Clear();

        }


        //convert the hashSet to array so max value can be extracted
        int repeatCount = totalRepeats.Count;
        int[] repeats = new int[repeatCount];
        totalRepeats.CopyTo(repeats);
        totalRepeats.Clear();

        //convert vertices matrix to vector array so rows can be deleted
        List<Vector<double>> vert = new List<Vector<double>>();
        for (int i = 0; i < numVertices; i++)
        {
            vert.Add(v.Row(i));
        }

        //delete all the repeated indices from vertices array
        for (int i = 0; i < repeatCount; i++)
        {
            int max = repeats.Max();
            vert.RemoveAt(max);
            repeats = repeats.Where(val => val != max).ToArray(); //removing max value from the array 'repeats'
        }

        int newNumVertices = numVertices - repeatCount;
        //Debug.Log(newNumVertices);

        Matrix<double> new_v = Matrix<double>.Build.Dense(newNumVertices, 3);

        for (int i = 0; i < (numVertices - repeatCount); i++)
        {
            new_v.SetRow(i, vert[i]);
        }

        //Debug.Log(new_v);


        Vector3[] newVertices = new Vector3[newNumVertices];

        for (int i = 0; i < newNumVertices; i++)  //fill in the vertices matrix
        {
            newVertices[i] = new Vector3((float)(new_v[i, 0]), (float)(new_v[i, 1]), (float)(new_v[i, 2]));
        }
        //Debug.Log(f);
        return f;



    }
}

