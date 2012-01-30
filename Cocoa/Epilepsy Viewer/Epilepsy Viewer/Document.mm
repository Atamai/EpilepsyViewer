//
//  Document.mm
//  Epilepsy Viewer
//
//  Created by Yves Starreveld on 12-01-23.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "Document.h"

#import "BasicVTKView.h"

// from EpilepsyViewer
#include "EpilepsyViewerData.h"
#include "EpilepsyViewerDisplay.h"
#include "EpilepsyViewerInteraction.h"

// from VTK
#import "vtkInteractorStyleSwitch.h"
#import "vtkCocoaRenderWindowInteractor.h"
#import "vtkConeSource.h"
#import "vtkCylinderSource.h"
#import "vtkPolyDataMapper.h"
#import "vtkSmartPointer.h"
#import "vtkDebugLeaks.h"

// from ToolCursor
#include "vtkToolCursorInteractorObserver.h"

@implementation Document

- (void)setupMainVTKView
{
    [mainVTKView initializeVTKSupport];
    
    // data, display and interaction managers
    EpilepsyViewerData dataManagerObject;
    EpilepsyViewerDisplay displayManagerObject;
    EpilepsyViewerInteraction interactionManagerObject;
    
    // get pointers to the objects (for stylistic reasons only)
    EpilepsyViewerData *dataManager = &dataManagerObject;
    EpilepsyViewerDisplay *displayManager = &displayManagerObject;
    EpilepsyViewerInteraction *interactionManager = &interactionManagerObject;

    NSString *shortPath = @"/Users/ystarrev/Desktop/AHN_MINC";
    char *buffer = new char[512];
    [shortPath getCString:buffer];

    dataManager->SetDataDirectory(buffer);

    if (!dataManager->LoadFromDataDirectory())
    {
        NSLog(@"Unable to read data, exiting.\n");
    }

    
    // go through the data processing steps
    dataManager->RegisterMRHeadToCTHead();
    dataManager->ExtractMRBrain();
    dataManager->ExtractCTElectrodes();
    
    // push the data to the display manager
    displayManager->SetData(dataManager);
        
    // for "basic" epilepsy viewer, use interactor observer
    vtkSmartPointer<vtkToolCursorInteractorObserver> observer =
    vtkSmartPointer<vtkToolCursorInteractorObserver>::New();
    observer->SetToolCursor(interactionManager->GetToolCursor());
    observer->SetInteractor([mainVTKView getInteractor]);
    observer->SetEnabled(1);
    
    displayManager->SetRenderWindow((vtkRenderWindow *)[mainVTKView getVTKRenderWindow]);
    interactionManager->BindInteraction(displayManager->GetMainRenderer());
    
     
/*    vtkSmartPointer<vtkConeSource> cone = vtkSmartPointer<vtkConeSource>::New();
    cone->SetHeight(3.0);
    cone->SetRadius(1.0);
    cone->SetResolution(100);
    vtkSmartPointer<vtkPolyDataMapper>  coneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    coneMapper->SetInput(cone->GetOutput());
    vtkSmartPointer<vtkActor>   coneActor = vtkSmartPointer<vtkActor>::New();
    coneActor->SetMapper(coneMapper);
    [mainVTKView getRenderer]->AddActor(coneActor);
    
    // Tell the system that the view needs to be redrawn
    [mainVTKView setNeedsDisplay:YES];
 */
}

- (id)init
{
    self = [super init];
    if (self) {
        // Add your subclass-specific initialization here.
        // If an error occurs here, return nil.
    }
    return self;
}

- (NSString *)windowNibName
{
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers, you should remove this method and override -makeWindowControllers instead.
    return @"Document";
}

- (void)windowControllerDidLoadNib:(NSWindowController *)aController
{
    [super windowControllerDidLoadNib:aController];
    // Add any code here that needs to be executed once the windowController has loaded the document's window.
    [self setupMainVTKView];

}

- (NSData *)dataOfType:(NSString *)typeName error:(NSError **)outError
{
    /*
     Insert code here to write your document to data of the specified type. If outError != NULL, ensure that you create and set an appropriate error when returning nil.
    You can also choose to override -fileWrapperOfType:error:, -writeToURL:ofType:error:, or -writeToURL:ofType:forSaveOperation:originalContentsURL:error: instead.
    */
    NSException *exception = [NSException exceptionWithName:@"UnimplementedMethod" reason:[NSString stringWithFormat:@"%@ is unimplemented", NSStringFromSelector(_cmd)] userInfo:nil];
    @throw exception;
    return nil;
}

- (BOOL)readFromData:(NSData *)data ofType:(NSString *)typeName error:(NSError **)outError
{
    /*
    Insert code here to read your document from the given data of the specified type. If outError != NULL, ensure that you create and set an appropriate error when returning NO.
    You can also choose to override -readFromFileWrapper:ofType:error: or -readFromURL:ofType:error: instead.
    If you override either of these, you should also override -isEntireFileLoaded to return NO if the contents are lazily loaded.
    */
    NSException *exception = [NSException exceptionWithName:@"UnimplementedMethod" reason:[NSString stringWithFormat:@"%@ is unimplemented", NSStringFromSelector(_cmd)] userInfo:nil];
    @throw exception;
    return YES;
}

+ (BOOL)autosavesInPlace
{
    return YES;
}

@end
